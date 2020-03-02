function eps_control_vector = control_vector(...
    t, nx, controlType, controlPressure, min_on, max_on,...
    lower_control_threshold, upper_control_threshold, Ii,...
    spatialControlLoc...
)
    switch controlType
        case 0 %No Control
            eps_control_vector = zeros(1,nx);         
        case 1 %Always-Everywhere
            eps_control_vector(1:nx) = controlPressure;                             
        case 2 %Time Dependent Control
            if(min_on < t && t < max_on)
                eps_control_vector(1:nx) = controlPressure;                                                            
            else
                eps_control_vector = zeros(1,nx);                       
            end
        case 3 %Periodic Control
            if(any(min_on < t & t < max_on))
                eps_control_vector(1:nx) = controlPressure;
            else
                eps_control_vector = zeros(1,nx);  
            end
            %use find(min_on < t & t < max_on) to determine application #
        case 4 %Density Dependent Control
            if(min_on < t && t < max_on)
                densityControlSwitch = zeros(1,nx);
                densityControlSwitch(Ii>lower_control_threshold&Ii<upper_control_threshold)=1;
                eps_control_vector = densityControlSwitch*controlPressure;
            else
                eps_control_vector = zeros(1,nx);     
            end
        case 5  %Spatial
            eps_control_vector = zeros(1,nx); 
            eps_control_vector(spatialControlLoc>0) = controlPressure;       
        otherwise %Implement No Control
            eps_control_vector = zeros(1,nx);       
    end            
end