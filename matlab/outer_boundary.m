function fb = outer_boundary(tr)
    
    map = dictionary();
    map("r") = "r";
    map("r") = [];
    for k=1:size(tr.ConnectivityList, 1)
        for i=1:2
        id_i = tr.ConnectivityList(k, i);
            for j=i+1:3
                id_j = tr.ConnectivityList(k, j);

                key = toString(min(id_i, id_j), max(id_i, id_j));

                value = toString(id_i, id_j);
                if (i==1 && j==3)
                    value = toString(id_j, id_i);
                end
                if (map.isKey(key))
                    map(key) = [];
                else
                    map(key) = value;
                end
            end
        end
    end
    
    ls = {};
    edges = map.values();
    for k = 1:size(edges, 1)
        edge =  fromString(edges(k));
        r1 = tr.Points(edge(1), 1);
        r2 = tr.Points(edge(2), 1);
        if (r1~=0 || r2~=0)
            ls{end+1} = edge;
        end
    end

    fb = zeros(size(ls, 1), 2);
    
    for k = 1:size(ls, 2)
        edge = ls{k};
        fb(k, 1) = edge(1);
        fb(k, 2) = edge(2);
    end

   
end

function s = toString(i,j)
s = string(i) + "," + string(j);
end
function edge = fromString(s)
    a = s.split(",");
    edge = [str2num(a(1)), str2num(a(2))];
end