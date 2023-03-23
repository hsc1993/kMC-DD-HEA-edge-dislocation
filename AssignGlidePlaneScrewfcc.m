function glideplane=AssignGlidePlaneScrewfcc(normal)

perfect_plane=[1 1 1; 1 1 -1;-1 1 1;1 -1 1]/norm([1 1 1]);
num_perf_planes=size(perfect_plane,1);
dot_first=abs(dot(perfect_plane(1,:),normal));
glideplane=perfect_plane(1,:);

for i=2:num_perf_planes
    check=abs(dot(perfect_plane(i,:),normal));
    if(check>dot_first)
        glideplane=perfect_plane(i,:);
        dot_first=check;
    end
end