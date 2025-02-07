function [threadXY, threadXYVel, hThread ,u] = ThreadUpdateBarrier(threadXY, threadXYVel, timeStep, hThread,u ,duC)


threadXY = threadXY + timeStep*threadXYVel;
u = u + timeStep*duC;


hThread.XData = [u(1,:) ,threadXY(1,:)];
hThread.YData = [u(2,:) ,threadXY(2,:)];

end