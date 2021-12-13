function MAX = decrepeak(y)

for i=1:length(y)
    if y(i)<0
        break;
    end
end
j =i;
t = 1;
while j<length(y)
    if y(j)>0
        start = j;
        while j<length(y)
            if y(j)<0
                last = j;
                MAX(t) = max(y(start:last-1));
                t = t+1;
                j = j+1;
                break;
            else
                j = j+1;
                continue
            end
        end
    else
        j = j+1;
    end
end
end
        
    