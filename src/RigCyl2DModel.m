classdef RigCyl2DModel < Fluid & RigidCylinder2D & SpringDamper2D
    %RigCylModel should compose the three class above, instead of being a
    %subclass of them
    
    properties
        MassRatio_x
        MassRatio_y
        NominalNaturalFreq_x
        NominalNaturalFreq_y
    end
    
    methods
        function obj = RigCyl2DModel(diameter, length, density_x, density_y, fluidSpeed, ...
                fluidDensity, springCons_x, springCons_y, damperCons_x, damperCons_y )
          
          obj = obj@RigidCylinder2D();  
          obj = obj@Fluid();
          obj = obj@SpringDamper2D();
          
          if(nargin ==10)   
                obj.Diameter = diameter;
                obj.Length = length;
                obj.Density_x = density_x;
                obj.Density_y = density_y;
                obj.FluidSpeed = fluidSpeed; %m/s
                obj.FluidDensity = fluidDensity; %kg/m^3
                obj.SpringCons_x = springCons_x;
                obj.DamperCons_x = damperCons_x;
                obj.SpringCons_y = springCons_y;
                obj.DamperCons_y = damperCons_y;
          else
                
            obj.Diameter = 0.1;
            obj.Length = 10;
            obj.Density_x = 1000;
            obj.Density_y = 1000;
            
            obj.FluidSpeed = 1; %m/s
            obj.FluidDensity = 1000; %kg/m^3
            
            obj.SpringCons_x = 100;
            obj.DamperCons_x = 100;
            obj.SpringCons_y = 100;
            obj.DamperCons_y = 100;
            
          end
          
          obj.Mass_x = obj.Density_x* pi* obj.Diameter^2 * obj.Length / 4;
          obj.Mass_y = obj.Density_y* pi* obj.Diameter^2 * obj.Length / 4;
          obj.MassRatio_x = obj.Mass_x/(obj.FluidDensity* pi* obj.Diameter^2 * obj.Length / 4);
          obj.MassRatio_y = obj.Mass_y/(obj.FluidDensity* pi* obj.Diameter^2 * obj.Length / 4);
           
          obj.NominalNaturalFreq_x = sqrt(obj.SpringCons_x/(obj.Mass_x * ...
               (1 + 1/obj.MassRatio_x)))/2/pi; 
          obj.NominalNaturalFreq_y = sqrt(obj.SpringCons_y/(obj.Mass_y * ...
               (1 + 1/obj.MassRatio_y)))/2/pi;            
        end
        
        function PrintRigCylModel2D(obj)
            disp(' ');
            disp('System Settings: ');
            disp(' ');
            disp(['Mass Ratio Inline = ', num2str(obj.MassRatio_x)]);
            disp(['Nominal Natural Frequence Inline = ', num2str(obj.NominalNaturalFreq_x) ' Hz']); 
            disp(['Nominal NonDimensional Natural Frequence Inline = ', num2str(obj.FluidSpeed/obj.Diameter/obj.NominalNaturalFreq_x)]); 
            disp(['Mass Ratio Crossflow = ', num2str(obj.MassRatio_y)]);
            disp(['Nominal Natural Frequence Crossflow = ', num2str(obj.NominalNaturalFreq_y) ' Hz']); 
            disp(['Nominal NonDimensional Natural Frequence Crossflow = ', num2str(obj.FluidSpeed/obj.Diameter/obj.NominalNaturalFreq_y)]); 
            disp(' ');
            disp('Cylinder Property:');
            disp(['Diameter = ', num2str(obj.Diameter), ' m']);
            disp(['Length = ', num2str(obj.Length), ' m']);
            disp(['Density Inline = ', num2str(obj.Density_x), ' kg/m^3']);
            disp(['Density Crossflow = ', num2str(obj.Density_y), ' kg/m^3']);
            disp(['Mass Inline = ', num2str(obj.Mass_x), ' kg']);
            disp(['Mass Crossflow = ', num2str(obj.Mass_y), ' kg']);
            disp(' ');
            disp('Fluid Property:');
            disp(['Fluid Speed = ', num2str(obj.FluidSpeed), ' m/s']);
            disp(['Fluid Density = ', num2str(obj.FluidDensity), ' kg/m^3']);
            disp(' ');
            disp('Spring-Damper Property:');
            disp(['Spring Constant Inline= ', num2str(obj.SpringCons_x)]);
            disp(['Damper Constant Inline= ', num2str(obj.DamperCons_x)]);
            disp(['Spring Constant Crossflow= ', num2str(obj.SpringCons_y)]);
            disp(['Damper Constant Crossflow= ', num2str(obj.DamperCons_y)]);
           
        end
       
%         function obj = set.FluidSpeed(obj, fluidSpeed)
%             obj.FluidSpeed = fluidSpeed;
%         end

        
    end
    
end

