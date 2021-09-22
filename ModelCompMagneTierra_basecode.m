classdef ModelCompMagneTierra_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel  matlab.ui.control.Label
        NmerodepartculasEditFieldLabel  matlab.ui.control.Label
        NmerodepartculasEditField       matlab.ui.control.NumericEditField
        ValorfinaldelrangoEditFieldLabel  matlab.ui.control.Label
        ValorfinaldelrangoEditField     matlab.ui.control.NumericEditField
        ValorinicialdelrangoEditField_2Label  matlab.ui.control.Label
        ValorinicialdelrangoEditField_2  matlab.ui.control.NumericEditField
        RangodevelocidadesparalaspartculasLabel  matlab.ui.control.Label
        ComenzarSimulacinButton         matlab.ui.control.Button
        ParticulasparaSimulacinLabel    matlab.ui.control.Label
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
        
        function [xFin, vFin]= eulerModificado2(app, q,m, xInic, vInic, delta, cx, cy ,cz, I, dominio)
            n = 30;
            % ------------------- aceleracion inicial con Lorentz
            [Bi,Bj,Bk] = campoMagnetico3D(app,xInic(1),xInic(2),xInic(3),cx,cy,cz,I,n,dominio);
            
            comp_i = vInic(2)*Bk - vInic(3)*Bj;
            comp_j = - vInic(1)*Bk + vInic(3)*Bi;
            comp_k = - vInic(2)*Bi + vInic(1)*Bj;
        
            aInic = (q/m)*[comp_i, comp_j, comp_k];
            
            % ------------------- EULER 1
            xFinEuler = xInic + delta * vInic;
            vFinEuler = vInic + delta * aInic;
            
            %   producto cruz v x B
            [Bi,Bj,Bk] = campoMagnetico3D(app,xFinEuler(1),xFinEuler(2),xFinEuler(3),cx,cy,cz,I,n,dominio);
            
            comp_i = vFinEuler(2)*Bk - vFinEuler(3)*Bj;
            comp_j = - vFinEuler(1)*Bk + vFinEuler(3)*Bi;
            comp_k = - vFinEuler(2)*Bi + vFinEuler(1)*Bj;
            
            comp_cruz = [comp_i, comp_j, comp_k];
            aFinEuler = (q/m)*comp_cruz;
            
            % --------------------- EULER MODIFICADO
            
            vFin = vInic + delta *((aInic + aFinEuler)/2);
            xFin = xInic + delta * ((vInic +vFinEuler)/2);
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
            
            % calcluar las integrales de las componentes de campo magnetico 
            
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
        
        function [p]= coordParti(app, q,m, xInic, vInic, delta, cx, cy ,cz, I, dominio)
            partx = [];
            party = [];
            partz = [];
            
            for i=1:delta:30
                [xFin, vFin]= eulerModificado2(app,q,m, xInic, vInic, delta, cx, cy ,cz, I, dominio);
                partx = [partx, xFin(1)];
                party = [party, xFin(2)];
                partz = [partz, xFin(3)];
                
                xInic = xFin;
                vInic = vFin;
            end
            p = [partx; party; partz];
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
    
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ComenzarSimulacinButton
        function ComenzarSimulacinButtonPushed(app, event)
            %% Programa Magnetico 
            % Entragable del Reto 
            % Profesor: Mario Iván Estrada Delgado
            % Profesor: Francisco Javier Delgado Cepeda
            
            % EQUIPO 1
            % Ana Patricia Islas Mainou (A01751676)
            % Pablo González de la Parra (A01745096)
            % Paula Sophia Santoyo Arteaga (A01745312)
            % Valeria Martínez Silva (A01752167)
            % Ricardo Antonio Cervantes Martínez (A01745912)
            
            n = 6;
            xmin = -200; xmax = 200; ymin = -200; ymax = 200; zmin = -200; zmax= 200;
            I = 0.001;
            
            %parametrizaciones
            syms s
            
            % esfera
            dominio = [-99,99];
            cx = 100*sqrt(1-(0.01*s)^2)*cos(s);
            cy = 100*sqrt(1-(0.01*s)^2)*sin(s);
            cz = 100*0.01*s;
            
            % genera mallado de coordenadas iniciales 
            [posx, posy, posz] = dominioRes3D(app,xmin, xmax, ymin, ymax, zmin, zmax, n);
            
            % genera campo magnetico 
            [Mx,My,Mz] = campoMagnetico3D(app,posx,posy,posz,cx,cy,cz,I,n,dominio);
            
            % posicion de las particulas
            delta = 0.1; 
           
            % se da por hecho que todas las particulas son electrones
            q = -1.6*10^-19; m = 9.11*10^-31;
            
            % genera particulas aleatorias
            nPuntos = app.NmerodepartculasEditField.Value; % numero de particulas
        
            % "y" y "z" cambian en aleaotrio, "x" se queda fijo}
            % rango sobre x
            a = xmin + 10; %xmin
            b = xmax - 10;  %xmax
        
            X = xmin +1;
            
            Y = randi([a,b],1,nPuntos); % vector de números aleatorios entre -50 y 50
            Z = randi([a,b],1,nPuntos); % vector de números aleatorios entre -50 y 50
        
            % las velocidades "vy" y "vz" se quedan en 0 y "vx" se aleatorio
            % rango para la velocidad
            a = app.ValorinicialdelrangoEditField_2.Value;
            b = app.ValorfinaldelrangoEditField.Value;
            
            vX = randi([a,b],1,nPuntos); % vector de números aleatorios entre -50 y 50
            vY = 0;
            vZ = 0;
        
            p_xyz = [];
            velInic = [];
        
            for i = 1:length(Y)
                p_xyz = [p_xyz; X, Y(i), Z(i)];
                velInic = [velInic; vX(i), vY, vZ];
            end
            
            [filas, col] = size(p_xyz);
            for i = 1:filas
                
                 p = coordParti(app,q,m, p_xyz(i,:), velInic(i,:), delta, cx, cy ,cz, I, dominio);
                 pos(:,:,i)= p;
            end
            
            % graficar
            
            xlim(app.UIAxes,[xmin,xmax])
            ylim(app.UIAxes,[ymin,ymax])
            zlim(app.UIAxes,[zmin,zmax])
            
            graficador3D(app,posx, posy, posz, Mx, My, Mz);
            hold(app.UIAxes ,"on")
            fplot3(app.UIAxes,cx,cy,cz,[dominio(1), dominio(2)],"b")
            
            hold(app.UIAxes ,"off")
            
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
            app.UIFigure.Color = [0.6235 0.8588 0.3373];
            app.UIFigure.Position = [100 100 590 611];
            app.UIFigure.Name = 'MATLAB App';

            % Create SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel
            app.SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel = uilabel(app.UIFigure);
            app.SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel.HorizontalAlignment = 'center';
            app.SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel.WordWrap = 'on';
            app.SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel.FontName = 'Arial';
            app.SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel.FontSize = 24;
            app.SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel.FontWeight = 'bold';
            app.SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel.Position = [53 527 485 73];
            app.SimulacindeunaTormentadeionesenelcampomagnticodelatierraLabel.Text = 'Simulación de una Tormenta de iones en el campo magnético de la tierra';

            % Create NmerodepartculasEditFieldLabel
            app.NmerodepartculasEditFieldLabel = uilabel(app.UIFigure);
            app.NmerodepartculasEditFieldLabel.HorizontalAlignment = 'right';
            app.NmerodepartculasEditFieldLabel.Position = [338 63 120 22];
            app.NmerodepartculasEditFieldLabel.Text = 'Número de partículas';

            % Create NmerodepartculasEditField
            app.NmerodepartculasEditField = uieditfield(app.UIFigure, 'numeric');
            app.NmerodepartculasEditField.Position = [473 63 55 22];
            app.NmerodepartculasEditField.Value = 1;

            % Create ValorfinaldelrangoEditFieldLabel
            app.ValorfinaldelrangoEditFieldLabel = uilabel(app.UIFigure);
            app.ValorfinaldelrangoEditFieldLabel.HorizontalAlignment = 'right';
            app.ValorfinaldelrangoEditFieldLabel.FontName = 'Arial';
            app.ValorfinaldelrangoEditFieldLabel.Position = [87 18 111 22];
            app.ValorfinaldelrangoEditFieldLabel.Text = 'Valor final del rango';

            % Create ValorfinaldelrangoEditField
            app.ValorfinaldelrangoEditField = uieditfield(app.UIFigure, 'numeric');
            app.ValorfinaldelrangoEditField.Position = [213 18 55 22];
            app.ValorfinaldelrangoEditField.Value = 45;

            % Create ValorinicialdelrangoEditField_2Label
            app.ValorinicialdelrangoEditField_2Label = uilabel(app.UIFigure);
            app.ValorinicialdelrangoEditField_2Label.HorizontalAlignment = 'right';
            app.ValorinicialdelrangoEditField_2Label.FontName = 'Arial';
            app.ValorinicialdelrangoEditField_2Label.Position = [79 47 119 22];
            app.ValorinicialdelrangoEditField_2Label.Text = 'Valor inicial del rango';

            % Create ValorinicialdelrangoEditField_2
            app.ValorinicialdelrangoEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.ValorinicialdelrangoEditField_2.Position = [213 47 55 22];
            app.ValorinicialdelrangoEditField_2.Value = 25;

            % Create RangodevelocidadesparalaspartculasLabel
            app.RangodevelocidadesparalaspartculasLabel = uilabel(app.UIFigure);
            app.RangodevelocidadesparalaspartculasLabel.HorizontalAlignment = 'center';
            app.RangodevelocidadesparalaspartculasLabel.WordWrap = 'on';
            app.RangodevelocidadesparalaspartculasLabel.FontName = 'Arial';
            app.RangodevelocidadesparalaspartculasLabel.FontSize = 16;
            app.RangodevelocidadesparalaspartculasLabel.FontWeight = 'bold';
            app.RangodevelocidadesparalaspartculasLabel.Position = [79 84 189 37];
            app.RangodevelocidadesparalaspartculasLabel.Text = 'Rango de velocidades para las partículas';

            % Create ComenzarSimulacinButton
            app.ComenzarSimulacinButton = uibutton(app.UIFigure, 'push');
            app.ComenzarSimulacinButton.ButtonPushedFcn = createCallbackFcn(app, @ComenzarSimulacinButtonPushed, true);
            app.ComenzarSimulacinButton.Position = [367 18 133 22];
            app.ComenzarSimulacinButton.Text = 'Comenzar Simulación';

            % Create ParticulasparaSimulacinLabel
            app.ParticulasparaSimulacinLabel = uilabel(app.UIFigure);
            app.ParticulasparaSimulacinLabel.HorizontalAlignment = 'center';
            app.ParticulasparaSimulacinLabel.WordWrap = 'on';
            app.ParticulasparaSimulacinLabel.FontName = 'Arial';
            app.ParticulasparaSimulacinLabel.FontSize = 16;
            app.ParticulasparaSimulacinLabel.FontWeight = 'bold';
            app.ParticulasparaSimulacinLabel.Position = [327 95 211 26];
            app.ParticulasparaSimulacinLabel.Text = 'Particulas para Simulación';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            app.UIAxes.Position = [65 128 462 389];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ModelCompMagneTierra_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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