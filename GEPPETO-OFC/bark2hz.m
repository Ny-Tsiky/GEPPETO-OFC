function res = bark2hz(f)
lomask = f<2;
himask = f>20.1;
f(lomask) = (f(lomask)-0.3)/0.85;
f(himask) = (f(himask)+4.422)/1.22;
res = 1960*(f+0.53)./(26.28-f);
end

% [m,n] = size(f);
% res = zeros(m,n);
% for j = 1:n
%     for i = 1:m
%         if (f(i,j) < 2)
%             f(i,j) = (f(i,j)-0.3)/0.85;
%         end
%         if (f(i,j) > 20.1)
%             f(i,j) = (f(i,j)+4.422)/1.22;
%         end
%     end
%     res(:,j) = 1960*((f(:,j)+0.53)./(26.28-f(:,j)));
% end
% end
