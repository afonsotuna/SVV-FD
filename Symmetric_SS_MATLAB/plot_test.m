t_lookup = 3635;
t_limit = 30; %seconds after

index = find(time == t_lookup);
t = 0:0.1:t_limit;
step_y = t_limit/0.1;
y = y_impulse(index:index+step_y,1);
clf()
plot(t,y)