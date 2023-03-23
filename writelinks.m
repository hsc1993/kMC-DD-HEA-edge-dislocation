function [] = writelinks(links)
%Writelinks takes the links file and writes it to a text file

ftemp = fopen('links.txt','W');

%Printing format
fmt = ["%f %f %f\n"];

for i = 1:length(links)
    fprintf(ftemp,fmt,links(i,1:3));
end

fclose(ftemp);

end

 