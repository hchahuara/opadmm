function colors = generate_contrast_colors(n)
    % Generate n distinct colors with good contrast
    hues = linspace(0, 1, n+1)'; % Generate hues, extra one to avoid wrap-around overlap
    hues(end) = []; % Remove the last one to keep exactly n hues
    saturations = repmat([1; 0.5], ceil(n/2), 1); % Alternate between full and half saturation
    saturations = saturations(1:n); % Ensure we only take n saturations
    values = repmat([1; 0.7], ceil(n/2), 1); % Alternate between full and reduced value
    values = values(1:n); % Ensure we only take n values
    
    hsv_colors = [hues, saturations, values];
    rgb_colors = hsv2rgb(hsv_colors);
    colors_hex = cellstr(rgb2hex(rgb_colors));
    
 colors = char(colors_hex); % Convert cell array to char matrix
end