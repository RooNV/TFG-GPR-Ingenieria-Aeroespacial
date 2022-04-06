function ans = Epoch2seconds(Epoch1, Epoch2)
    dv = datevec(Epoch1);
    dt = datetime(dv,'Format','dd-MM-yyyy HH:mm:ss.SSSSSS');
    dv2 = datevec(Epoch2);
    dt2 = datetime(dv2,'Format','dd-MM-yyyy HH:mm:ss.SSSSSS');
    ans = seconds(dt-dt2);
end