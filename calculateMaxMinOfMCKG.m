minM = min(abs(Parameter.Matrix.mass(find(Parameter.Matrix.mass~=0))));
maxM = max(abs(Parameter.Matrix.mass(find(Parameter.Matrix.mass~=0))));
minC = min(abs(Parameter.Matrix.damping(find(Parameter.Matrix.damping~=0))));
maxC = max(abs(Parameter.Matrix.damping(find(Parameter.Matrix.damping~=0))));
minK = min(abs(Parameter.Matrix.stiffness(find(Parameter.Matrix.stiffness~=0))));
maxK = max(abs(Parameter.Matrix.stiffness(find(Parameter.Matrix.stiffness~=0))));
minG = min(abs(Parameter.Matrix.gyroscopic(find(Parameter.Matrix.gyroscopic~=0))));
maxG = max(abs(Parameter.Matrix.gyroscopic(find(Parameter.Matrix.gyroscopic~=0))));
values = [maxM, maxC, maxK, maxG;...
          minM, minC, minK, minG];
values = full(values);
ratio = values(1,:) ./ values(2,:);
