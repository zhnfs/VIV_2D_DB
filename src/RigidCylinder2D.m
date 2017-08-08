classdef RigidCylinder2D
    %RIGIDCYLINDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Diameter
        Length
        Density_x
        Density_y
        Mass_x
        Mass_y
    end
    
    methods
        function obj = RigidCylinder2D(diameter, length, Mass_x, Mass_y)
            if(nargin == 3)   
                obj.Diameter = diameter;
                obj.Length = length;
                obj.Mass_x = Mass_x;
                obj.Mass_y = Mass_y;
                obj.Density_x = obj.Mass_x/(pi*obj.Diameter^2*obj.Length/4);
                obj.Density_y = obj.Mass_y/(pi*obj.Diameter^2*obj.Length/4);
            else
                obj.Diameter = 0.1;
                obj.Length = 10;
                obj.Density_x = 1000;
                obj.Density_y = 1000;
                obj.Mass_x = obj.Density_x*pi*obj.Diameter^2*obj.Length/4;
                obj.Mass_y = obj.Density_y*pi*obj.Diameter^2*obj.Length/4;
            end
            
        end


%         end
    end
    
end

