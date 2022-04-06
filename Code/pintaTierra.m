function pintaTierra()
    % Pinta la Tierra
    [x,y,z] = sphere(50);
    load topo;
    props.AmbientStrength = 0.1;
    props.DiffuseStrength = 1;
    props.SpecularColorReflectance = .3; %.5;
    props.SpecularExponent = 20;
    props.SpecularStrength = 1;
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo;
    surface(6700000*x,6700000*y,6700000*z,props); %radio
    light('position',[-1 0 1]); % -1 0 1
    light('position',[1.5 -0.5 -0.5], 'color', [.6 .2 .2]); % -1.5 .5 -.5
    view([-90 0]) % az = -90 y el = 0
    axis equal off
    hold on
end