function make_it_nicer
    axisProperties=gca; % get the pointer to the axis properties
    axisProperties.FontSize=15; % set the font size
    axisProperties.Title.Interpreter='latex'; % set the text interpreter
    axisProperties.Subtitle.Interpreter='latex'; % set the text interpreter
    axisProperties.YLabel.Interpreter='latex';
    axisProperties.XLabel.Interpreter='latex';
end