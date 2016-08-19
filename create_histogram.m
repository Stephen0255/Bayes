function [ ] = create_histogram( model, parameters, several_scheme, BinWidth )
% Create the histograms for the severals parameters

if several_scheme
    figure()
    histogram(parameters.scheme)
    title('Distribution of the schemes');
    xlabel('Downwind, Center, Upwind');

    if strcmp(model,'first_order_homogenous')
        a_upwind=parameters.coefficients_upwind(1,:);
        a_downwind=parameters.coefficients_downwind(1,:);
        a_center=parameters.coefficients_center(1,:);
    
        figure()
        histogram(a_upwind,'BinWidth', BinWidth)
        title('Upwind');
        xlabel('a');
    
        figure()
        histogram(a_center,'BinWidth', BinWidth)
        title('Center');
        xlabel('a');
    
        figure()
        histogram(a_downwind,'BinWidth', BinWidth)
        title('Downwind');
        xlabel('a');
        
    elseif strcmp(model,'first_order')    
        a_upwind=parameters.coefficients_upwind(1,:);
        a_downwind=parameters.coefficients_downwind(1,:);
        a_center=parameters.coefficients_center(1,:);
        b_upwind=parameters.coefficients_upwind(2,:);
        b_downwind=parameters.coefficients_downwind(2,:);
        b_center=parameters.coefficients_center(2,:);    
    
        figure()
        histogram(a_upwind,'BinWidth', BinWidth)
        title('Upwind');
        xlabel('a');
        
        figure()
        histogram(b_upwind,'BinWidth', BinWidth)
        title('Upwind');
        xlabel('b');   
    
        figure()
        histogram(a_center,'BinWidth', BinWidth)
        title('Center');
        xlabel('a');
        
        figure()
        histogram(b_center,'BinWidth', BinWidth)
        title('Center');
        xlabel('b');
    
        figure()
        histogram(a_downwind,'BinWidth', BinWidth)
        title('Downwind');
        xlabel('a');
        
        figure()
        histogram(b_downwind,'BinWidth', BinWidth)
        title('Downwind');
        xlabel('b');    
    end
    
elseif ~several_scheme

    if parameters.scheme == 1
        str = 'Upwind';
    elseif parameters.scheme == 0
        str = 'Center';
    elseif  parameters.scheme == -1
        str = 'Downwind';
    end
    
    if strcmp(model,'first_order_homogenous')
        figure()
        histogram(parameters.coefficients(1,:),'BinWidth', BinWidth)
        title(str);
        xlabel('a');
        
    elseif strcmp(model,'first_order')
        figure()
        histogram(parameters.coefficients(1,:),'BinWidth', BinWidth)
        title(str);
        xlabel('a');
        
        figure()
        histogram(parameters.coefficients(2,:),'BinWidth', BinWidth)
        title(str);
        xlabel('b');
    end
end


end

