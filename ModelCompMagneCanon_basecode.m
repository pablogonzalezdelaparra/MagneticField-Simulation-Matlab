classdef ModelCompMagneCanon_basecode < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        SimulacindelcandeelectronesLabel  matlab.ui.control.Label
        SolenoidearribaEditFieldLabel   matlab.ui.control.Label
        SolenoidearribaEditField        matlab.ui.control.NumericEditField
        IntensidaddecorrientedelossolenoidesLabel  matlab.ui.control.Label
        SolenoideabajoEditFieldLabel    matlab.ui.control.Label
        SolenoideabajoEditField         matlab.ui.control.NumericEditField
        SolenoidederechaEditFieldLabel  matlab.ui.control.Label
        SolenoidederechaEditField       matlab.ui.control.NumericEditField
        SolenoideizquierdaEditFieldLabel  matlab.ui.control.Label
        SolenoideizquierdaEditField     matlab.ui.control.NumericEditField
        InformacindelaparticulaLabel    matlab.ui.control.Label
        VelocidadInicialsobreXEditFieldLabel  matlab.ui.control.Label
        VelocidadInicialsobreXEditField  matlab.ui.control.NumericEditField
        ComenzarSimulacinButton         matlab.ui.control.Button
        UIAxes                          matlab.ui.control.UIAxes
    end

    
    methods (Access = private)
        
        function area= integraTrapezoidal(app, funcion, n, a, b)
            
            paso = (b-a)/n;
            
            % coordenadas x
            x = a:paso:b;
            
            % coordenades y f(x)
            fx = funcion(x);
            
            % usa patron 1 2 2 ... 2 2 1
            fx = [fx(1), fx(2:length(fx)-1) .*2 , fx(length(fx))] ;
            
            % calcula area bajo la curva
            area =  sum(fx) * (b-a)/(2*n);
        end
        
        function graficador3D(app, x, y, z, Fx, Fy, Fz)
            % Generera las magnitudes
            [a,b,c] = size(Fx);
            mag= zeros(a,a,a);
            for i=1:1:a
                  for j = 1:1:b
                      
                      for k = 1:1:c
                            mag(i,j,k) = sqrt(Fx(i,j,k)^2 + Fy(i,j,k)^2 + Fz(i,j,k)^2);
                      end
                  end
            end
        
            
            % re escala los vectores unitarios a 0.25
            unit_x = (Fx ./mag).*40;
            unit_y = (Fy ./mag).*40;
            unit_z = (Fz ./mag).*40;
            
            % obtener opacidad
            %op = 1 - rescale(mag,0.3,1);
            mag = (((mag)/(max(max(max(mag)))))*0.7)+ 0.3;
            op = 1 - mag;
            
            op(isnan(op))= 0;
            
                hold(app.UIAxes, "on")
            for i=1:1:a
                  for j = 1:1:b
                      for k = 1:1:c
                         
                      quiver3(app.UIAxes, x(i,j), y(i,j), z(k),unit_x(i,j,k), unit_y(i,j,k),unit_z(i,j,k), 0, 'Color', ...
                          [op(i,j,k), op(i,j,k), op(i,j,k)], "AutoScale", 0, "MaxHeadSize",3);
        
                      end
                  end
            end
                 
                 hold(app.UIAxes, "on")
                 view(app.UIAxes, 40,40)
                end
        
        function [x, y, z] = dominioRes3D(app, x1, x2, y1, y2,z1, z2, n)
            % para la matrix x
            x = [];
            for i = 1:n
            x = [x; linspace(x1,x2,n)];
            end
            
            % para la matrix y 
            y =[];
            for i = 1:n
            y = [y; linspace(y1,y2,n)];
            end
            y = transpose(y);
            
            % para el vector en z
            z = linspace(z1,z2,n);
        end
        
        function [Mx,My, Mz] = campoMagnetico3D(app,x,y,z,cx,cy,cz,I,n,dominio)
            [r,q] = size(x);
    
            Mx = zeros(r,r,r);
            My = zeros(r,r,r);
            Mz = zeros(r,r,r);
            
            % constantes necesarias
            mo = 4*pi*10^-7; %permiabilidad magnetica
            
            syms s
            % obtener dl = dcds/ds no ponemos el ds de arriba aqui
            dlx = diff(cx); %componente i
            dly = diff(cy); %componente j
            dlz = diff(cz); %componente k
            
            % calcular las integrales de las componentes de campo magnetico 
            
            for i = 1:r
                for j = 1:r
                    for k = 1:r
                        
                        % funciones que se van a integrar
                       [intX, intY, intZ] = fIntegral(app,x(i,j),y (i,j),z(k),cx,cy,cz,dlx, dly, dlz);
                       try
                            Mx(i,j,k) = integraTrapezoidal(app,intX, n, dominio(1), dominio(2));
                            My(i,j,k) = integraTrapezoidal(app,intY, n, dominio(1), dominio(2));
                            Mz(i,j,k) = integraTrapezoidal(app,intZ, n, dominio(1), dominio(2));
                       catch
                           
                            Mx(i,j,k) = 0;
                            My(i,j,k) = 0;
                            Mz(i,j,k) = 0;
                       end
                        
                    end
                end
            end
            
            Mx = Mx .*((mo*I)/(4*pi));
            My = My .*((mo*I)/(4*pi));
            Mz = Mz .*((mo*I)/(4*pi));
        end
        
        function [intX, intY, intZ] =fIntegral(app,x,y,z,cx,cy,cz,dlx, dly, dlz)
            syms s
    
            rx = x-cx; %componente i
            ry = y-cy; %componente j
            rz = z-cz; %componente k
            
            % producto cruz (pC) dlxr
            pC_x = dly*rz - dlz*ry; %componente i
            pC_y = - dlx*rz + dlz*rx; %componente j
            pC_z = dlx*ry - dly*rx; %componente k
            
            
            intX(s) = (pC_x/(rx^2 + ry^2 + rz^2) ^(3/2));
            intY(s) = (pC_y/(rx^2 + ry^2 + rz^2) ^(3/2));
            intZ(s) = (pC_z/(rx^2 + ry^2 + rz^2) ^(3/2));
        end
        
        function trayectoria (app,x,y,z)
            
            h= animatedline(app.UIAxes,"Marker", "o", "Linestyle", "-","Color", "r", "MarkerSize", 5, "MarkerFaceColor", "r");
            h2 = animatedline(app.UIAxes,"Linestyle", "-", "Color","r", "MarkerFaceColor", "r");
            
            for i=1:length(x)
                clearpoints(h)
                addpoints(h, x(i), y(i), z(i));
                addpoints(h2, x(i), y(i), z(i));
                
                drawnow
                view(app.UIAxes,40,40)
            end
        end
        
        function particulas(app,pos)
            
            tam = size(pos); %%num de particulas en 3era componente
            tam_obj = size(size(pos));
            
            if tam_obj(2) == 2
                % cuando hay una particula
                trayectoria(app,pos(1,:), pos(2,:), pos(3,:))
                
            else
                % cuando hay mas de 1 particula
                npart = tam(3);
                for i=1:npart
                    f(i)= animatedline(app.UIAxes,"Marker", "o", "Linestyle", "-", "Color","r", "MarkerSize", 5, "MarkerFaceColor", "r");
                    g(i)= animatedline(app.UIAxes,"Linestyle", "-", "Color","r", "MarkerSize", 5, "MarkerFaceColor", "r");
                end
                nstep = tam(2);
        
                for i=1:nstep
                    for j= 1:npart
                        clearpoints(f(j));
                        addpoints(f(j), pos(1,i,j), pos(2,i,j),pos(3,i,j));
                        addpoints(g(j), pos(1,i,j), pos(2,i,j),pos(3,i,j));
                        
                        drawnow
                        view(app.UIAxes,40,40)
                    end
                end
            end
        end
        
        function [xFin, vFin]= eulerModificado2canon(app,q,m, xInic, vInic, delta, cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3,cx4, cy4,cz4,I1,I2,I3,I4, dominio)
            n = 30;
            % ------------------- aceleracion inicial con Lorentz
            [Bi1,Bj1,Bk1] = campoMagnetico3D(app,xInic(1),xInic(2),xInic(3),cx1,cy1,cz1,I1,n,dominio);
            [Bi2,Bj2,Bk2] = campoMagnetico3D(app,xInic(1),xInic(2),xInic(3),cx2,cy2,cz2,I2,n,dominio);
            [Bi3,Bj3,Bk3] = campoMagnetico3D(app,xInic(1),xInic(2),xInic(3),cx3,cy3,cz3,I3,n,dominio);
            [Bi4,Bj4,Bk4] = campoMagnetico3D(app,xInic(1),xInic(2),xInic(3),cx4,cy4,cz4,I4,n,dominio);
            
            Bi = Bi1 + Bi2 + Bi3 + Bi4;
            Bj = Bj1 + Bj2 + Bj3 + Bj4;
            Bk = Bk1 + Bk2 + Bk3 + Bk4;
            
            comp_i = vInic(2)*Bk - vInic(3)*Bj;
            comp_j = - vInic(1)*Bk + vInic(3)*Bi;
            comp_k = - vInic(2)*Bi + vInic(1)*Bj;
        
            aInic = (q/m)*[comp_i, comp_j, comp_k];
            
            % ------------------- EULER 1
            xFinEuler = xInic + delta * vInic;
            vFinEuler = vInic + delta * aInic;
            
            %   producto cruz v x B
            [Bi1,Bj1,Bk1] = campoMagnetico3D(app,xFinEuler(1),xFinEuler(2),xFinEuler(3),cx1,cy1,cz1,I1,n,dominio);
            [Bi2,Bj2,Bk2] = campoMagnetico3D(app,xFinEuler(1),xFinEuler(2),xFinEuler(3),cx2,cy2,cz2,I2,n,dominio);
            [Bi3,Bj3,Bk3] = campoMagnetico3D(app,xFinEuler(1),xFinEuler(2),xFinEuler(3),cx3,cy3,cz3,I3,n,dominio);
            [Bi4,Bj4,Bk4] = campoMagnetico3D(app,xFinEuler(1),xFinEuler(2),xFinEuler(3),cx4,cy4,cz4,I4,n,dominio);
            
            Bi = Bi1 + Bi2 + Bi3 + Bi4;
            Bj = Bj1 + Bj2 + Bj3 + Bj4;
            Bk = Bk1 + Bk2 + Bk3 + Bk4;
            
            comp_i = vFinEuler(2)*Bk - vFinEuler(3)*Bj;
            comp_j = - vFinEuler(1)*Bk + vFinEuler(3)*Bi;
            comp_k = - vFinEuler(2)*Bi + vFinEuler(1)*Bj;
            
            comp_cruz = [comp_i, comp_j, comp_k];
            aFinEuler = (q/m)*comp_cruz;
            
            % --------------------- EULER MODIFICADO
            
            vFin = vInic + delta *((aInic + aFinEuler)/2);
            xFin = xInic + delta * ((vInic +vFinEuler)/2);
    
        end
        
        function [p]= coordPartiCanon(app,q,m, xInic, vInic, delta, cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3,cx4, cy4,cz4, I1,I2,I3,I4, dominio)

            partx = [];
            party = [];
            partz = [];
            
            for i=1:delta:2
                [xFin, vFin]= eulerModificado2canon(app,q,m, xInic, vInic, delta, cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3,cx4, cy4,cz4,I1,I2,I3,I4, dominio);
                partx = [partx, xFin(1)];
                party = [party, xFin(2)];
                partz = [partz, xFin(3)];
                
                xInic = xFin;
                vInic = vFin;
            end
            p = [partx; party; partz];
                end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            hold(app.UIAxes ,"on")
            
            % resorte
            syms s
            
            dominio = [-15,15];
            
            % solenoide 1 arriba
            cx1 = 100*cos(s);
            cy1 = 100*sin(s);
            cz1 = 100*(1.7 + 0.05*s);
            
            % solenoide 2 abajo
            cx2 = 100*cos(s);
            cy2 = 100*sin(s);
            cz2 = 100*(-1.7 + 0.05*s);
            
            % solenoide 3 izquierda
            cx3 = 100*cos(s);
            cz3 = 100*(sin(s));
            cy3 = 100*(1.7 + 0.05*s);
            
            % solenoide 4 derecha
            cx4 = 100*cos(s);
            cz4 = 100*(sin(s));
            cy4 = 100*(-1.7 + 0.05*s);
            
            fplot3(app.UIAxes,cx1,cy1,cz1,[dominio(1), dominio(2)],"b");
            fplot3(app.UIAxes,cx2,cy2,cz2,[dominio(1), dominio(2)],"b");
            fplot3(app.UIAxes,cx3,cy3,cz3,[dominio(1), dominio(2)],"g");
            fplot3(app.UIAxes,cx4,cy4,cz4,[dominio(1), dominio(2)],"g");
            
            view(app.UIAxes,40,40)
            hold(app.UIAxes ,"off")
        end

        % Button pushed function: ComenzarSimulacinButton
        function ComenzarSimulacinButtonPushed(app, event)
            cla(app.UIAxes) 
            
            n = 6;
            xmin = -150; xmax = 150; ymin = -300; ymax =300; zmin = -350; zmax= 350;
            
            % resorte
            syms s
            
            dominio = [-15,15];
            
            % solenoide 1 arriba
            cx1 = 100*cos(s);
            cy1 = 100*sin(s);
            cz1 = 100*(1.7 + 0.05*s);
            I1 = app.SolenoidearribaEditField.Value;
            
            % solenoide 2 abajo
            cx2 = 100*cos(s);
            cy2 = 100*sin(s);
            cz2 = 100*(-1.7 + 0.05*s);
            I2 = app.SolenoideabajoEditField.Value;
            
            % solenoide 3 izquierda
            cx3 = 100*cos(s);
            cz3 = 100*(sin(s));
            cy3 = 100*(1.7 + 0.05*s);
            I3 = app.SolenoideizquierdaEditField.Value;
            
            % solenoide 4 derecha
            cx4 = 100*cos(s);
            cz4 = 100*(sin(s));
            cy4 = 100*(-1.7 + 0.05*s);
            I4 = app.SolenoidederechaEditField.Value;
            
            % genera mallado de coordenadas iniciales 
            [posx, posy, posz] = dominioRes3D(app,xmin, xmax, ymin, ymax, zmin, zmax, n);
            
            % campos de cada solenoide
            [Mx1,My1,Mz1] = campoMagnetico3D(app,posx,posy,posz,cx1,cy1,cz1,I1,n,dominio);
            [Mx2,My2,Mz2] = campoMagnetico3D(app,posx,posy,posz,cx2,cy2,cz2,I2,n,dominio);
            [Mx3,My3,Mz3] = campoMagnetico3D(app,posx,posy,posz,cx3,cy3,cz3,I3,n,dominio);
            [Mx4,My4,Mz4] = campoMagnetico3D(app,posx,posy,posz,cx4,cy4,cz4,I4,n,dominio);
            
            % campo total 
            Mx = Mx1 + Mx2 + Mx3 + Mx4;
            My = My1 + My2 + My3 + My4;
            Mz = Mz1 + Mz2 + Mz3 + Mz4;
            
            % posicion de las particulas
            delta = 0.05; 
           
            % se da por hecho que todas las particulas son electrones
            q = -1.6*10^-19; m = 9.11*10^-31;
            
            X = -150;
            Y = 0;
            Z = 0;
            
            vX = app.VelocidadInicialsobreXEditField.Value;
            vY = 0;
            vZ = 0;
            
            p_xyz = [];
            velInic = [];
            
            for i = 1:length(Y)
                p_xyz = [p_xyz; X, Y(i), Z(i)];
                velInic = [velInic; vX(i), vY, vZ];
            end
            
            [filas, col]= size(p_xyz);
            for i = 1:filas
                
                 p = coordPartiCanon(app,q,m, p_xyz(i,:), velInic(i,:), delta, cx1,cy1,cz1,cx2,cy2,cz2,cx3,cy3,cz3,cx4, cy4,cz4, I1,I2,I3,I4, dominio);
                 pos(:,:,i)= p;
            end
            
            %    GRAFICAR

            xmin = xmin -50; xmax = xmax +50;
            ymin = ymin -50; ymax = ymax +50;
            zmin = zmin -50; zmax = zmax +50;
            
            xlim(app.UIAxes ,[xmin,xmax])
            ylim(app.UIAxes ,[ymin,ymax])
            zlim(app.UIAxes ,[zmin,zmax])
            
            hold(app.UIAxes ,"on")
            graficador3D(app,posx, posy, posz, Mx, My, Mz);
            hold(app.UIAxes ,"on")
            
            fplot3(app.UIAxes,cx1,cy1,cz1,[dominio(1), dominio(2)],"b");
            fplot3(app.UIAxes,cx2,cy2,cz2,[dominio(1), dominio(2)],"b");
            fplot3(app.UIAxes,cx3,cy3,cz3,[dominio(1), dominio(2)],"g");
            fplot3(app.UIAxes,cx4,cy4,cz4,[dominio(1), dominio(2)],"g");
            
            view(app.UIAxes,60,30)
             
            hold(app.UIAxes ,"on")
            particulas(app,pos);
               
            hold(app.UIAxes ,"off")
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.3608 0.6275 0.8784];
            app.UIFigure.Position = [100 100 640 467];
            app.UIFigure.Name = 'MATLAB App';

            % Create SimulacindelcandeelectronesLabel
            app.SimulacindelcandeelectronesLabel = uilabel(app.UIFigure);
            app.SimulacindelcandeelectronesLabel.HorizontalAlignment = 'center';
            app.SimulacindelcandeelectronesLabel.WordWrap = 'on';
            app.SimulacindelcandeelectronesLabel.FontName = 'Arial';
            app.SimulacindelcandeelectronesLabel.FontSize = 24;
            app.SimulacindelcandeelectronesLabel.FontWeight = 'bold';
            app.SimulacindelcandeelectronesLabel.Position = [78 406 485 41];
            app.SimulacindelcandeelectronesLabel.Text = 'Simulación del cañón de electrones';

            % Create SolenoidearribaEditFieldLabel
            app.SolenoidearribaEditFieldLabel = uilabel(app.UIFigure);
            app.SolenoidearribaEditFieldLabel.HorizontalAlignment = 'right';
            app.SolenoidearribaEditFieldLabel.FontName = 'Arial';
            app.SolenoidearribaEditFieldLabel.Position = [458 307 93 22];
            app.SolenoidearribaEditFieldLabel.Text = 'Solenoide arriba';

            % Create SolenoidearribaEditField
            app.SolenoidearribaEditField = uieditfield(app.UIFigure, 'numeric');
            app.SolenoidearribaEditField.Position = [566 307 46 22];
            app.SolenoidearribaEditField.Value = 0.001;

            % Create IntensidaddecorrientedelossolenoidesLabel
            app.IntensidaddecorrientedelossolenoidesLabel = uilabel(app.UIFigure);
            app.IntensidaddecorrientedelossolenoidesLabel.HorizontalAlignment = 'center';
            app.IntensidaddecorrientedelossolenoidesLabel.WordWrap = 'on';
            app.IntensidaddecorrientedelossolenoidesLabel.FontName = 'Arial';
            app.IntensidaddecorrientedelossolenoidesLabel.FontSize = 16;
            app.IntensidaddecorrientedelossolenoidesLabel.FontWeight = 'bold';
            app.IntensidaddecorrientedelossolenoidesLabel.Position = [438 342 189 38];
            app.IntensidaddecorrientedelossolenoidesLabel.Text = 'Intensidad de corriente de los solenoides';

            % Create SolenoideabajoEditFieldLabel
            app.SolenoideabajoEditFieldLabel = uilabel(app.UIFigure);
            app.SolenoideabajoEditFieldLabel.HorizontalAlignment = 'right';
            app.SolenoideabajoEditFieldLabel.FontName = 'Arial';
            app.SolenoideabajoEditFieldLabel.Position = [460 277 91 22];
            app.SolenoideabajoEditFieldLabel.Text = 'Solenoide abajo';

            % Create SolenoideabajoEditField
            app.SolenoideabajoEditField = uieditfield(app.UIFigure, 'numeric');
            app.SolenoideabajoEditField.Position = [566 277 46 22];
            app.SolenoideabajoEditField.Value = 0.001;

            % Create SolenoidederechaEditFieldLabel
            app.SolenoidederechaEditFieldLabel = uilabel(app.UIFigure);
            app.SolenoidederechaEditFieldLabel.HorizontalAlignment = 'right';
            app.SolenoidederechaEditFieldLabel.FontName = 'Arial';
            app.SolenoidederechaEditFieldLabel.Position = [446 246 105 22];
            app.SolenoidederechaEditFieldLabel.Text = 'Solenoide derecha';

            % Create SolenoidederechaEditField
            app.SolenoidederechaEditField = uieditfield(app.UIFigure, 'numeric');
            app.SolenoidederechaEditField.Position = [566 246 46 22];
            app.SolenoidederechaEditField.Value = 0.001;

            % Create SolenoideizquierdaEditFieldLabel
            app.SolenoideizquierdaEditFieldLabel = uilabel(app.UIFigure);
            app.SolenoideizquierdaEditFieldLabel.HorizontalAlignment = 'right';
            app.SolenoideizquierdaEditFieldLabel.FontName = 'Arial';
            app.SolenoideizquierdaEditFieldLabel.Position = [440 216 111 22];
            app.SolenoideizquierdaEditFieldLabel.Text = 'Solenoide izquierda';

            % Create SolenoideizquierdaEditField
            app.SolenoideizquierdaEditField = uieditfield(app.UIFigure, 'numeric');
            app.SolenoideizquierdaEditField.Position = [566 216 46 22];
            app.SolenoideizquierdaEditField.Value = 0.001;

            % Create InformacindelaparticulaLabel
            app.InformacindelaparticulaLabel = uilabel(app.UIFigure);
            app.InformacindelaparticulaLabel.HorizontalAlignment = 'center';
            app.InformacindelaparticulaLabel.WordWrap = 'on';
            app.InformacindelaparticulaLabel.FontName = 'Arial';
            app.InformacindelaparticulaLabel.FontSize = 16;
            app.InformacindelaparticulaLabel.FontWeight = 'bold';
            app.InformacindelaparticulaLabel.Position = [441 140 189 38];
            app.InformacindelaparticulaLabel.Text = 'Información de la particula';

            % Create VelocidadInicialsobreXEditFieldLabel
            app.VelocidadInicialsobreXEditFieldLabel = uilabel(app.UIFigure);
            app.VelocidadInicialsobreXEditFieldLabel.HorizontalAlignment = 'right';
            app.VelocidadInicialsobreXEditFieldLabel.WordWrap = 'on';
            app.VelocidadInicialsobreXEditFieldLabel.FontName = 'Arial';
            app.VelocidadInicialsobreXEditFieldLabel.Position = [447 93 111 31];
            app.VelocidadInicialsobreXEditFieldLabel.Text = 'Velocidad Inicial sobre X';

            % Create VelocidadInicialsobreXEditField
            app.VelocidadInicialsobreXEditField = uieditfield(app.UIFigure, 'numeric');
            app.VelocidadInicialsobreXEditField.Position = [573 97 46 22];
            app.VelocidadInicialsobreXEditField.Value = 600;

            % Create ComenzarSimulacinButton
            app.ComenzarSimulacinButton = uibutton(app.UIFigure, 'push');
            app.ComenzarSimulacinButton.ButtonPushedFcn = createCallbackFcn(app, @ComenzarSimulacinButtonPushed, true);
            app.ComenzarSimulacinButton.Position = [469 46 133 22];
            app.ComenzarSimulacinButton.Text = 'Comenzar Simulación';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [18 24 415 383];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ModelCompMagneCanon_basecode

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end