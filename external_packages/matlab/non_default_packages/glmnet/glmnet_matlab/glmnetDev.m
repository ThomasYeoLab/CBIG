function dev = glmnetDev(object)

dev = (1-object.dev) * object.nulldev;

end