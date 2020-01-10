set(gca, 'defaultUicontrolFontName', 'Ubuntu')
set(gca, 'defaultUitableFontName', 'Ubuntu')
set(gca, 'defaultAxesFontName', 'Ubuntu')
set(gca, 'defaultTextFontName', 'Ubuntu')
set(gca, 'defaultUipanelFontName', 'Ubuntu')
fig = gcf;
set(findall(fig, '-property', 'FontSize'), 'FontSize', 16)