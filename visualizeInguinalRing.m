function [penetrationPoins,inguinalAbdominalWallPolygons,hInguinalRing] = visualizeInguinalRing(varargin)
%visualizeInguinalRing2 Plots a patch object using the geometry of the
%inguinal ring

penetrationPoins =1;

parser = inputParser;
parser.addParameter('NeedleRadius', 0.000355);
parser.addParameter('FractionEdgeToRing', 0.5);
parse(parser, varargin{:})

needleRadius = parser.Results.NeedleRadius;
fractionEdgeToRing = parser.Results.FractionEdgeToRing;

linPar = linspace(0,1,11);

%%% Define distances
innerRad  = 0.009/2*1000*2; 
outerRad  = 0.006*1000*2;
curveRad  = 0.0044*1000*2; 
sideLen   = 0.01*1000*2;
botTotLen = 0.06*1000*2;
topBotOuterLen = fractionEdgeToRing*(botTotLen/2 - (curveRad + outerRad) - needleRadius); % [m]
topInnerLen    = (1-fractionEdgeToRing)*(botTotLen/2 - (curveRad + outerRad)) - needleRadius; % [m]
botInnerLen    = botTotLen - 2*topBotOuterLen - 2*needleRadius - 2*needleRadius; % [m]

%%% Define initial vertices
% Left rectangle
LeftRectangleTopLeftCorner  = [-(outerRad + curveRad + topInnerLen + 2*needleRadius + topBotOuterLen)-0.5;-curveRad];
% Right rectangle
RightRectangleTopLeftCorner = [ (outerRad + curveRad + topInnerLen + 2*needleRadius)+0.5;-curveRad];
% RightRectangleTopLeftCorner = [ (outerRad + curveRad + topInnerLen + 2*needleRadius)+0;-curveRad];

% Center Piece
CenterPieceTopLeftCorner = [-(outerRad + curveRad + topInnerLen );-curveRad];

%%% Define the vertex edge segments
% Left Rectangle%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LeftRectangleTopLine   = LeftRectangleTopLeftCorner *[1,1] + (topBotOuterLen-0.001)*[ 1; 0]*[0,1];
LeftRectangleRightSide = LeftRectangleTopLine(:,end)   *[1,1] + sideLen       *[ 0;-1]*[0,1];
LeftRectangleBotLine   = LeftRectangleRightSide(:,end)*[1,1] + (topBotOuterLen-0.001)*[-1; 0]*[0,1];
LeftRectangleLeftSide  = LeftRectangleBotLine(:,end)  *[1,1] + sideLen       *[ 0; 1]*[0,1];
LeftRectangle = [LeftRectangleRightSide(:,1:end-1),LeftRectangleBotLine(:,1:end-1),LeftRectangleLeftSide(:,1:end-1),LeftRectangleTopLine(:,1:end-1)]; %2X4 double

% curving edges
pos = [LeftRectangleLeftSide(1,1:end-1) LeftRectangleLeftSide(2,1:end-1) LeftRectangleRightSide(1,1:end-1)-LeftRectangleLeftSide(1,1:end-1) LeftRectangleRightSide(2,1:end-1)-LeftRectangleLeftSide(2,1:end-1)]; % [x y w h]
cur = [1 1]*0.25; % curvature

npts = 20;
npts = round(npts/4)*4 + 1;
th = linspace(0,2*pi,npts);
thm = reshape(th(2:end),[],4);
thm = [th(1) thm(end,1:3); thm];

rr = pos(3:4)/2; %w/2 h/2
er = cur.*rr; % radii of the elliptical arcs 
x = er(1)*cos(thm) + pos(1) + rr(1) + (rr(1)-er(1))*[1 -1 -1 1];
y = er(2)*sin(thm) + pos(2) + rr(2) + (rr(2)-er(2))*[1 1 -1 -1];

x = [x(:); x(1)];
y = [y(:); y(1)];
LeftRectangle = [x'; y'];


% Right Rectangle%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RightRectangleTopLine   = RightRectangleTopLeftCorner *[1,1] + topBotOuterLen*[ 1; 0]*[0,1];
RightRectangleRightSide = RightRectangleTopLine(:,end)   *[1,1] + sideLen       *[ 0;-1]*[0,1];
RightRectangleBotLine   = RightRectangleRightSide(:,end)*[1,1] + topBotOuterLen*[-1; 0]*[0,1];
RightRectangleLeftSide  = RightRectangleBotLine(:,end)  *[1,1] + sideLen       *[ 0; 1]*[0,1];
RightRectangle = [RightRectangleRightSide(:,1:end-1),RightRectangleBotLine(:,1:end-1),RightRectangleLeftSide(:,1:end-1),RightRectangleTopLine(:,1:end-1)];

% curving edges
pos = [RightRectangleLeftSide(1,1:end-1) RightRectangleLeftSide(2,1:end-1) RightRectangleRightSide(1,1:end-1)-RightRectangleLeftSide(1,1:end-1) RightRectangleRightSide(2,1:end-1)-RightRectangleLeftSide(2,1:end-1)]; % [x y w h]

rr = pos(3:4)/2; %w/2 h/2
er = cur.*rr; % radii of the elliptical arcs 
x = er(1)*cos(thm) + pos(1) + rr(1) + (rr(1)-er(1))*[1 -1 -1 1];
y = er(2)*sin(thm) + pos(2) + rr(2) + (rr(2)-er(2))*[1 1 -1 -1];

x = [x(:); x(1)];
y = [y(:); y(1)];
RightRectangle = [x'; y'];

% Center Piece%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CenterPieceTopLeftLine = CenterPieceTopLeftCorner      *[1,1] + topInnerLen     *[1;0]*[0,1];
CenterPieceLeftCurve   = CenterPieceTopLeftLine(:,end) + curveRad*[0;1] + curveRad*[cosd(270+90*linPar);sind(270+90*linPar)];
CenterPieceTopCurve   = CenterPieceLeftCurve(:,end) + outerRad*[1;0] + outerRad*[cosd(180-180*linPar);sind(180-180*linPar)];
CenterPieceRightCurve   = CenterPieceTopCurve(:,end) + curveRad*[1;0] + curveRad*[cosd(180+90*linPar);sind(180+90*linPar)];
CenterPieceTopRightLine = CenterPieceRightCurve(:,end) *[1,1] + topInnerLen     *[1;0]*[0,1];
CenterPieceRightLine = CenterPieceTopRightLine(:,end) *[1,1] + sideLen     *[0;-1]*[0,1];
CenterPieceBottomLine = CenterPieceRightLine(:,end) *[1,1] + botInnerLen     *[-1;0]*[0,1];
CenterPieceLeftLine = CenterPieceBottomLine(:,end) *[1,1] + sideLen     *[0;1]*[0,1];
% CenterPiece = [CenterPieceTopLeftLine(:,1:end-1),CenterPieceLeftCurve(:,1:end-1),CenterPieceTopCurve(:,1:end-1),CenterPieceRightCurve(:,1:end-1),CenterPieceTopRightLine(:,1:end-1),CenterPieceRightLine(:,1:end-1),CenterPieceBottomLine(:,1:end-1),CenterPieceLeftLine(:,1:end-1)];

% curving edges
pos = [CenterPieceLeftLine(1,1:end-1) CenterPieceLeftLine(2,1:end-1) CenterPieceRightLine(1,1:end-1)-CenterPieceLeftLine(1,1:end-1) CenterPieceRightLine(2,1:end-1)-CenterPieceLeftLine(2,1:end-1)]; % [x y w h]

rr = pos(3:4)/2; %w/2 h/2
er = cur.*rr; % radii of the elliptical arcs 
x = er(1)*cos(thm) + pos(1) + rr(1) + (rr(1)-er(1))*[1 -1 -1 1];
y = er(2)*sin(thm) + pos(2) + rr(2) + (rr(2)-er(2))*[1 1 -1 -1];

x = [x(:); x(1)];
y = [y(:); y(1)];
CenterPieceTopRightcurve = [x(1:6)';y(1:6)']; %from bottom right to top left
CenterPieceTopLeftcurve = [x(7:12)';y(7:12)'];  %from top right to bottom left
CenterPieceBottomLeftcurve = [x(13:18)';y(13:18)'];
CenterPieceBottomRightcurve = [x(19:25)';y(19:25)'];

CenterPiece = [flip(CenterPieceTopLeftcurve,2) , CenterPieceLeftCurve(:,2:end-1),CenterPieceTopCurve(:,1:end-1),CenterPieceRightCurve(:,1:end-1), flip(CenterPieceTopRightcurve,2), flip(CenterPieceBottomRightcurve(:,1:end-1),2) ,flip(CenterPieceBottomLeftcurve(:,1:end-1),2) ];

% Inner ring hole
InnerRing = innerRad*[cos(-linPar*2*pi);sin(-linPar*2*pi)];
% Convex hull
inguinalAbdominalWall = [LeftRectangleTopLine(:,1:end-1),CenterPieceLeftCurve(:,1:end-1),CenterPieceTopCurve(:,1:end-1),CenterPieceRightCurve(:,1:end-1),CenterPieceTopRightLine(:,1:end-1),RightRectangleRightSide(:,1:end-1),RightRectangleBotLine(:,1:end-1),LeftRectangleLeftSide(:,1:end-1)];
%%% Get Polygons
inguinalAbdominalWallWithNeedleHoles = [LeftRectangle,nan(2,1),CenterPiece,nan(2,1),RightRectangle];
inguinalAbdominalWallPolygons = {LeftRectangle.';CenterPiece.';RightRectangle.'};
polyInguinalRingInner = polyshape([inguinalAbdominalWallWithNeedleHoles,nan(2,1),InnerRing(:,1:end-1)].');
polyInguinalRingOuter = polyshape([inguinalAbdominalWall,nan(2,1),InnerRing(:,1:end-1)].');
%% Penetration Points
% penetrationPoins = [...
%     (CenterPieceBottomLine(:,1)+RightRectangleBotLine(:,2)).'/2;
%     (CenterPieceTopRightLine(:,2)+RightRectangleTopLine(:,1)).'/2;
%     (CenterPieceTopLeftLine(:,1)+LeftRectangleTopLine(:,2)).'/2;
%     (CenterPieceBottomLine(:,2)+LeftRectangleBotLine(:,1)).'/2;
%     ];

%% Visualize
% clf
%draws the pink thing
hInguinalRingInner = plot(polyInguinalRingInner,'LineStyle','none','FaceColor',[1 0.5 0.5],'FaceAlpha',0.5);
hold on
hInguinalRingConvHull = plot(polyInguinalRingOuter,'LineWidth',2,'FaceColor',[1 0.5 0.5],'EdgeColor',0.5*[1 0.5 0.5],'FaceAlpha',0.5);
hInguinalRing = [hInguinalRingInner;hInguinalRingConvHull];
% axis equal, grid on, grid minor
hInguinalRing =1;
end
