# Simple and efficient unzip, since apparently there is no Base.unzip yet
# NB. This actually returns a lazy iterator, so almost no upfront cost
unzip(a) = (getfield.(a, x) for x in fieldnames(eltype(a)))
