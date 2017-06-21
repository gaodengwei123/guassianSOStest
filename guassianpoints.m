function path = guassianpoints()
global Gaussiangao
path1 = [0         0
        4.07    0.1685
        23.12   24.79
    48.31   42.91
    50.0000   40.0000];

path2 = [0         0
    4.07    0.1685
    25.5236   14.4921
    46.5314   37.4165
    50.0000   40.0000];

path3 = [0         0
    1.777     14.38
    5.7     34.15
    29.46       47.41
    50.0000   40.0000];
path4 = [0         0
    1.777     14.38
    5.7     34.15
    19.38   37.41
    30.79   32.44
    50.0000   40.0000];
switch Gaussiangao
    case 1
        path = path1;
    case 2
        path = path2;
    case 3
        path = path3;
    case 4
        path = path4;
end



end



