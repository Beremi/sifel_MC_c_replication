% ************************************************************************
%
% transformation - computes nodal values from element values
%
% ************************************************************************

function node_value = transformation(elem_value)

  global elem nv nt

  node_value = zeros(1,nv);
  counter = zeros(1,nv);

  for jt = 1:nt
    ABC = elem(:,jt);
    node_value(ABC) = node_value(ABC) + elem_value(jt);
    counter(ABC) = counter(ABC) + 1;
  end

  active = counter > 0;
  node_value(active) = node_value(active)./counter(active);

end
