function [fitresult, gof] = create_fit_exp2(x, y)
%CREATEFIT1(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 3' fit:
%      X Input: x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 15-Apr-2024 11:02:23


%% Fit: 'untitled fit 3'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'exp2' );
excludedPoints = excludedata( xData, yData, 'Indices', [13 14 15 16 17 18 19 20 21 22 23 24 25 26 29 30 34 36 37 38 39 40 41 42 43 44 45 46 47 49 51 52 53 54 55 59 60 64 65 66 67 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 121 124 125 126 129 154 155 157 158 162 164 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 183 184 212 288 289 317 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 339 341 343 344 346 347 420 422 425 426 427 428] );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-10.0075542121138 0.481616965670719 11.5567522505614 -2.3386591803666];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 3' );
h = plot( fitresult, xData, yData, excludedPoints );
legend( h, 'y vs. x', 'Excluded y vs. x', 'untitled fit 3', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
grid on

