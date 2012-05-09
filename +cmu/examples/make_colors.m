
fid = fopen('color.database');
gid = fopen('colors.m','w');

fprintf(gid,'function c = colors(colorname)\n');
fprintf(gid,'% colors(''amber''\n');
fprintf(gid,'% returns RGB tuples adapted from\n');
fprintf(gid,'% http://en.wikipedia.org/wiki/List_of_colors\n');

fprintf(gid,'color_database = containers.Map();\n');

color_database = containers.Map();
while 1
    line = fgetl(fid);  
    if ~ischar(line), break, end
    if strcmp(line(1),'%'), continue, end % comment line
        
    t = regexp(line,'\t','split');
    [color,hex,r,g,b,hue,satur,light,satur2,value] = t{:};
    
    sr = strtrim(r);
    R = str2num(sr(1:end-1))/100;
    sg = strtrim(g);
    G = str2num(sg(1:end-1))/100;
    sb = strtrim(b);
    B = str2num(sb(1:end-1))/100;
    RGB = [R G B];
    %color_database(lower(strtrim(color))) = RGB;
    c = lower(strtrim(color));
    c = strrep(c,'''','');
    fprintf(gid,'color_database(''%s'')=[%f %f %f];\n',c,R,G,B)
end

fprintf(gid,'if strcmp(colorname,''list'')\n')
fprintf(gid,'k = color_database.keys;\n')
fprintf(gid,'k''\n')
fprintf(gid,'return; end\n')
fprintf(gid,'c = color_database(lower(colorname));');

fclose(fid)
fclose(gid)