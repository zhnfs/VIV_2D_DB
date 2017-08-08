classdef SpringDamper2D
    %SPRINGDASH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SpringCons_x = 100;
        SpringCons_y = 100;
        DamperCons_x = 100;
        DamperCons_y = 100;
    end
    
    methods
        function obj = SpringDamper2D(springCons_x, springCons_y, damperCons_x, damperCons_y)
          if( nargin == 2)
            obj.SpringCons_x = springCons_x;
            obj.DamperCons_x = damperCons_x;
            obj.SpringCons_y = springCons_y;
            obj.DamperCons_y = damperCons_y;
          else
            obj.SpringCons_x = 200;
            obj.DamperCons_x = 0;
            obj.SpringCons_y = 200;
            obj.DamperCons_y = 0;
          end
        end
    end
    
end

