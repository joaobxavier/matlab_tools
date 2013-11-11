function cmap = colormapOtu(otuNames)

otuOrderInColormap = {... 
'Enterobacteriaceae'
'Lactobacillus'
'Akkermansia'
'Barnesiella'
'Mollicutes'
'Clostridium_difficile'
'Blautia'
'Parasutterella'
'Turicibacter'
'Enterococcus'
'Bacteria'
'Allobaculum'
'Coprobacillus'
'Lachnospiraceae'
'Ruminococcaceae'
'Enterococcus'
'Ruminococcaceae'
'Oscillibacter'
'Clostridiales'
'other'};

cmapOtu = [...
254 147 33%     'Enterobacteriaceae'
0 128 255%     'Lactobacillus'
255 255 255%     'Akkermansia'
1 127 2%     'Barnesiella'
255 255 1%     'Mollicutes'
255 0 0%     'Clostridium_difficile'
205 102 255%     'Blautia'
255 204 102%     'Parasutterella'
65 0 128%     'Turicibacter'
34 60 251%     'Enterococcus'
0 0 0%     'Bacteria'
100 0 100%     'Allobaculum'
255 102 102%     'Coprobacillus'
1 255 255%     'Lachnospiraceae'
128 127 0%     'Ruminococcaceae'
0 0 254%     'Enterococcus'
125 125 125%     'Ruminococcaceae'
128 64 0%     'Oscillibacter'
125 0 0%     'Clostridiales'
50 50 50%     'other'
]/255;

if length(otuOrderInColormap) ~= size(cmapOtu, 1)
    error('number of otus is %d but colors is %d',...
        length(otuOrderInColormap), size(cmapOtu, 1));
end

cmap = zeros(length(otuNames), 3)+50/255;
for i = 1:length(otuOrderInColormap)
    otuIndeces = find(strcmp(otuOrderInColormap{i}, otuNames));
    if ~isempty(otuIndeces)
       for j = 1:length(otuIndeces)
           cmap(otuIndeces(j), :) = cmapOtu(i, :);
       end
    end
end