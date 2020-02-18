ks.wykreskorelacji = function(var1, var2, labvar1, labvar2, title, yx = T, metoda = 'pearson', gdzie_legenda = "topleft") {
  var1 = as.numeric(as.character(var1))
  var2 = as.numeric(as.character(var2))
  plot(var1, var2, pch = 19, cex=0.5, xlab=labvar1, ylab=labvar2, main=title)
  if (yx == T) { abline(0,1, col='gray') }
  abline(fit <- lm(var2 ~ var1), col='black', lty = 2)
  temp = cor.test(var1, var2, method = metoda)
  if (metoda=='pearson') {
    if (temp$p.value < 0.0001) {
      legend(gdzie_legenda, bty="n", legend=paste("r =",round(cor(var1,var2,use="complete.obs"),2),", p < 0.0001\nadj. R² =", format(summary(fit)$adj.r.squared, digits=4)))
    } else {
      legend(gdzie_legenda, bty="n", legend=paste("r =",round(cor(var1,var2,use="complete.obs"),2),", p =",round(temp$p.value, 4),"\nadj. R² =", format(summary(fit)$adj.r.squared, digits=4))) } }
  if (metoda=="spearman")
  { if (temp$p.value < 0.0001) {
    legend(gdzie_legenda, bty="n", legend=paste("rho =",round(cor(var1,var2,use="complete.obs"),2),", p < 0.0001"))
  } else {
    legend(gdzie_legenda, bty="n", legend=paste("rho =",round(cor(var1,var2,use="complete.obs"),2),", p =",round(temp$p.value, 4))) } } 
  
}

ks.correct01 = function(x, eng = T) {
  if (eng == T) {
    y = as.factor(ifelse(x==1, "1 Yes", "0 No")) }
  else { y = as.factor(ifelse(x==1, "1 Tak", "0 Nie")) }
  y
}

ks.correctvar = function(data, nvar, func) {
  newdata = do.call(func, list(data[, as.numeric(nvar)]))
  return(newdata)
}

ks.excelnumericdate = function(x)
{
  excel_numeric_to_date(as.numeric(as.character(x)), date_system = "modern") 
}

ks.plot_with_save = function(templot, tempname, height = 8.27, width = 11.69)
{
  #tiff(paste0("Plots/",tempname,".tiff"), res = 300, width = width, height = height)
  #templot
  #dev.off()
  #ggsave(paste0("Plots/",tempname,".tiff"), plot = templot, width = width, height = height, dpi = 300)
  templot
}