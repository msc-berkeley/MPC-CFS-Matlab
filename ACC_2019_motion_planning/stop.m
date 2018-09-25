rosinit('http://10.42.0.1:11311')

%%
a=1;
  while a<20
    pub = rospublisher('/cmd_vel','geometry_msgs/Twist');
    %pause(1)
    cmv = rosmessage('geometry_msgs/Twist');
    cmv.Linear.X = -0.0;
    cmv.Angular.Z = 0.0;
    send(pub,cmv);
    %pause(0.1)
    a=a+1;
  end
  
   
rosshutdown