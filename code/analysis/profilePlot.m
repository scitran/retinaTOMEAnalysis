function h=profilePlot(xVal, yVal, yMean, plotxlabel, plotylabel, plotTitle,plotFlag)

if plotFlag
    h=figure;
    plot(xVal,yVal,'-r');
    hold on
    plot(xVal,yMean,'-k','LineWidth',4);
    xlabel(plotxlabel);
    ylabel(plotylabel);
    title(plotTitle)
end