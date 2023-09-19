library(sf) # sudo apt install libudunits2-dev
library(sp)
library(automap)
library(maps)

dirname_src = utils::getSrcDirectory(function(){})[1]
if (is.na(dirname_src)) {
    dirname_src = "."
}
setwd(dirname_src)
source("./data_loading_and_merging.r")

LOAD_AND_DEFINE_REGION = function(fname_region_kml="../res/landscape.kml", fname_phenotype="../res/phenotype_data.csv", herbicide="Glyphosate", n_point_samples=1e4) {
    # fname_region_kml="../res/landscape.kml"
    # fname_phenotype="../res/phenotype_data.csv"
    # herbicide="Glyphosate"
    # n_point_samples=5000
    landscape = sf::st_read(fname_region_kml)
    phenotypes = LOAD_PHENOTYPES(fname_phenotype=fname_phenotype, phenotype_names=herbicide)
    phenotypes = phenotypes[!is.na(phenotypes$Longitude), ]
    phenotypes = phenotypes[!is.na(phenotypes$Latitude), ]
    x_limit = range(phenotypes$Longitude)
    y_limit = range(phenotypes$Latitude)
    x = phenotypes$Longitude
    y = phenotypes$Latitude
    z = eval(parse(text=paste0("phenotypes$", herbicide=herbicide)))
    idx = !is.na(z)
    x = x[idx]
    y = y[idx]
    z = z[idx]
    df_data = data.frame(x, y, z)
    ### Transform data.frame into and sp::SpatialPointsDataFrame
    df_data = aggregate(z ~ x + y, data=df_data, FUN=mean)
    sp::coordinates(df_data) = ~ x + y
    ### Prepare the new_data to fit the kriging model built using the above training data
    df_region = as.data.frame(as.matrix(landscape$geometry[[1]]))
    colnames(df_region) = c("x", "y", "z")
    sp::coordinates(df_region) = ~ x + y ### Transform data.frame into and sp::SpatialPointsDataFrame
    ### Prepare the new_data to fit the kriging model built using the above training data
    A = sp::Polygon(df_region) ### Transform into a polygon
    ### Sample regular intervals across the polygon to approximate its area and shape
    set.seed(42069) ### setting a randomisation seed for repeatability
    B = sp::spsample(A, n=n_point_samples, type="regular")
    return(list(df_data=df_data,
                new_data=B)
          )
}

KRIGING_AND_LINEAR_MODELS = function(REGION) {
    # REGION = LOAD_AND_DEFINE_REGION()
    vec_models = c("Ste", "Mat", "Sph", "Exp", "Gau")
    K1 = tryCatch(autoKrige.cv(z ~ 1, model=vec_models[1], REGION$df_data), error=function(e){autoKrige.cv(z ~ 1, REGION$df_data)})
    K2 = tryCatch(autoKrige.cv(z ~ 1, model=vec_models[2], REGION$df_data), error=function(e){autoKrige.cv(z ~ 1, REGION$df_data)})
    K3 = tryCatch(autoKrige.cv(z ~ 1, model=vec_models[3], REGION$df_data), error=function(e){autoKrige.cv(z ~ 1, REGION$df_data)})
    K4 = tryCatch(autoKrige.cv(z ~ 1, model=vec_models[4], REGION$df_data), error=function(e){autoKrige.cv(z ~ 1, REGION$df_data)})
    K5 = tryCatch(autoKrige.cv(z ~ 1, model=vec_models[5], REGION$df_data), error=function(e){autoKrige.cv(z ~ 1, REGION$df_data)})
    K_compare = compare.cv(K1, K2, K3, K4, K5)
    rmse = unlist(K_compare[rownames(K_compare)=="RMSE", ])
    idx = tail(which(rmse == min(rmse)), 1) # pick the last one, i.e. simpler models
    # mae = unlist(K_compare[rownames(K_compare)=="MAE", ])
    # idx = which(mae == min(mae))[1]
    model = vec_models[idx]
    K = tryCatch(automap::autoKrige(z ~ 1, REGION$df_data, model=model, new_data=REGION$new_data),
                    error=function(e) {automap::autoKrige(z ~ 1, REGION$df_data, new_data=REGION$new_data)}
    )
    ### Prepare the kriging output for plotting
    if (model=="Exp") {
        model = "Exponential"
    } else if (model=="Sph") {
        model = "Spherical"
    } else if (model=="Mat") {
        model = "Matern covariance"
    } else if (model=="Ste") {
        model = "Stein's Matern Parameterisation"
    } else {
        model = "Gaussian"
    }
    P = cbind(K$krige_output@coords, K$krige_output@data)[, 1:3]
    ### Set predicted resistances less than 0 to 0, and greater than 100 to 100
    P[P[,3] <   0, 3] = 0
    P[P[,3] > 100, 3] = 100
    colnames(P) = c("x", "y", "z") ### coordinates from the sampled points inside the paddock polygon and the krigng-predicted values, z
    ### Centre coordinates so that we can better account for the interaction effect
    x_centred = scale(REGION$df_data$x, center=TRUE, scale=FALSE)
    y_centred = scale(REGION$df_data$y, center=TRUE, scale=FALSE)
    z = REGION$df_data$z
    ### Model the herbicide gradient across the landscape
    mod = lm(z ~ x_centred + y_centred + (x_centred:y_centred))
    coef_summary = as.data.frame(summary(mod)$coefficients)
    sig = c()
    for (i in 1:nrow(coef_summary)){
        if (coef_summary[[4]][i] < 0.001) {
            sig = c(sig, "(***)")
        } else if (coef_summary[[4]][i] < 0.01) {
            sig = c(sig, "(**)")
        } else if (coef_summary[[4]][i] < 0.05) {
            sig = c(sig, "(*)")
        } else {
            sig = c(sig, "(ns)")
        }
    }
    coef_summary$sig_labels = paste(c("Intercept", "Longitude", "Latitude", "Interaction"), sig, sep="")
    return(list(df_kriged=P,
                model=model,
                coef_summary=coef_summary)
          )
}

PLOT_DISTRIBUTION_MAP = function(REGION, MODELS, fname_map_svg="../res/Glyphosate_resistance_distribution_SE_Australia.svg", label="Resistance (%)", 
    colours = list(rdylgn=rev(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#D9EF8B","#A6D96A","#66BD63","#1A9850","#006837")),
                   spectral=rev(c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2")))[[1]], 
    n_colours=101, plot_points=TRUE, plot_krig=TRUE, rescale_krig=FALSE, hist_not_range=TRUE) {

    # REGION = LOAD_AND_DEFINE_REGION()
    # MODELS = KRIGING_AND_LINEAR_MODELS(REGION)
    # fname_map_svg = "../res/Glyphosate_resistance_distribution_SE_Australia.svg"
    # n_colours = 101
    x = REGION$df_data$x
    y = REGION$df_data$y
    z = REGION$df_data$z
    x_limit = range(x)
    y_limit = range(y)
    vec_colours = colorRampPalette(colours)(n_colours)
    # width  = 6 * (diff(x_limit) / diff(y_limit))
    # height = 6 * 1
    width = 6 * (diff(x_limit) / diff(y_limit))
    height = 13 * (diff(y_limit) / diff(x_limit))
    svg(fname_map_svg, width=width, height=height)
    ### Plot the map
    plot(0, xlim=x_limit, ylim=y_limit, asp=1, type="n", xlab="Longitude", ylab="Latitude",
        main=paste0(gsub(".svg", "", gsub("_", " ", basename(fname_map_svg))),
                    "\n",
                    paste(tail(MODELS$coef_summary$sig_labels, 3), collapse=" | ")))
    grid()
    outline = maps::map("world", plot=FALSE)
    xrange = range(outline$x, na.rm=TRUE)
    yrange = range(outline$y, na.rm=TRUE)
    xbox = xrange + c(-2, 2)
    ybox = yrange + c(-2, 2)
    ### draw the outline of the map and color the water blue
    polypath(c(outline$x, NA, c(xbox, rev(xbox))),
            c(outline$y, NA, rep(ybox, each=2)),
            col="light blue", rule="evenodd")
    ### Plot kriging-predictions
    minimum = min(MODELS$df_kriged[, 3], na.rm=TRUE)
    maximum = max(MODELS$df_kriged[, 3], na.rm=TRUE)
    if (plot_krig==TRUE) {
        if (rescale_krig==TRUE) {
            if ((maximum - minimum) > 1e-7) {
                MODELS$df_kriged[, 3] = 100 * (MODELS$df_kriged[, 3] - minimum) / (maximum - minimum)
            } else {
                MODELS$df_kriged[, 3] = 50 ### set to mid colours if there is no variation in the kriged resistance levels
            }
            
        }
        for (k in 1:nrow(MODELS$df_kriged)){
            # k = 1
            idx = ceiling(MODELS$df_kriged[k, 3]) + 1
            points(MODELS$df_kriged$x[k], MODELS$df_kriged$y[k], pch=15, col=vec_colours[idx])
        }
        ### Append the best-fit Kriging model
        text(x=x_limit[1], y=y_limit[1], paste0("Kriging model:\n", MODELS$model), pos=4)
    }
    ### Plot populations and their resistance levels
    if (plot_points==TRUE) {
        for (i in 1:length(x)) {
            idx = ceiling(z[i]) + 1
            points(x[i], y[i], col="gray", bg=vec_colours[idx], pch=21)
        }
    }
    ### Heatmap legend
    nclass=10
    if (hist_not_range==TRUE) {
        par(fig=c(0.77,0.97,0.3,0.5), new=TRUE)
        par(mar=c(0,1,1,1))
        plot(0, ylab= "", xlab="", xaxt="n", yaxt="n", type="n")
        par(new=TRUE)
        df_table = data.frame(bins=seq(0, 100, length=11), counts=rep(0, times=11))
        for (i in 2:nrow(df_table)) {
            # i = 2
            df_table[i, 2] = sum((z >= df_table[i-1, 1]) & (z < df_table[i, 1]))
        }
        barplot(df_table$counts[2:nrow(df_table)], col=vec_colours[seq(1, n_colours, length=nclass)], bord=FALSE, las=1)
        par(fig=c(0.77,0.97,0.27,0.29), new=TRUE)
        par(mar=c(0,1,0,1))
        h = hist(seq(0,100, length=10), ylab= "", xlab="", xaxt="n", yaxt="n", las=1, main="", nclass=nclass, 
                    col=vec_colours[seq(1, n_colours, length=10)],
                    bord=FALSE)
        xrange = round(seq(h$breaks[1], h$breaks[length(h$breaks)], len=5), 2)
        labels = xrange
        axis(side=1, at=xrange, labels=labels, padj=-1)
        mtext(label, side=1, padj=2.5)
    } else {
        par(fig=c(0.77,0.97,0.40,0.50), new=TRUE)
        par(mar=c(1,1,1,1))
        h = hist(seq(0,100, length=10), ylab= "", xlab="", xaxt="n", yaxt="n", las=1, main="", nclass=nclass, 
                    col=vec_colours[seq(1, n_colours, length=10)],
                    bord=FALSE)
        xrange = round(seq(h$breaks[1], h$breaks[length(h$breaks)], len=5))
        labels = round(seq(minimum, maximum, len=5))
        axis(side=1, at=xrange, labels=labels, padj=-1)
        mtext(label, side=1, padj=2.5)
    }
    dev.off()
    return(fname_map_svg)
}

SAVE_SHAPEFILE = function(MODELS, fname_shapefile="../res/Glyphosate_resistance_gradient.shp", herbicide_name="Gly_resist") {
    ### Find the distances between sampled points
    dx = mean(diff(sort(unique(MODELS$df_kriged$x))))
    dy = mean(diff(sort(unique(MODELS$df_kriged$y))))
    ### Create qadrilaterals (squares in this case since we sampled in regular intervals from the A, i.e. the landscape)
    Q = lapply(1:nrow(MODELS$df_kriged), function(i){
        M = matrix(c(MODELS$df_kriged$x[i],    MODELS$df_kriged$y[i],
                    MODELS$df_kriged$x[i],    MODELS$df_kriged$y[i]+dy,
                    MODELS$df_kriged$x[i]+dx, MODELS$df_kriged$y[i]+dy,
                    MODELS$df_kriged$x[i]+dx, MODELS$df_kriged$y[i],
                    MODELS$df_kriged$x[i],    MODELS$df_kriged$y[i])
                , ncol =2, byrow = T
        )
        st_polygon(list(M))
    })
    eval(parse(text=paste0("SF_OBJECT = st_sf(", herbicide_name, "= MODELS$df_kriged$z, st_sfc(Q))")))
    sf::st_write(SF_OBJECT, fname_shapefile, append=FALSE) ### save into shapefile but the herbicide names will be trucated because of the inherent restrictions in the file format
    return(fname_shapefile)
}

HAVERSINE_DISTANCE_IN_KM = function(lon1, lat1, lon2, lat2, sphere_radius_km=6371) {
    ### Formula derved from: https://www.movable-type.co.uk/scripts/latlong.html
    y1 = lat1 * pi/180
    y2 = lat2 * pi/180
    x1 = lon1 * pi/180
    x2 = lon2 * pi/180
    dy = y2 - y1
    dx = x2 - x1
    a = sin(dy/2)^2 + (cos(y1) * cos(y2) * (sin(dx/2)^2))
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    d = sphere_radius_km * c
    return(d)
}

##################
### UNIT TESTS ###
##################
kriging_and_maps = function() {
    test_that(
        "LOAD_AND_DEFINE_REGION", {
            print("LOAD_AND_DEFINE_REGION:")
            REGION = LOAD_AND_DEFINE_REGION()
            expect_equal(dim(as.data.frame(REGION$df_data)), c(110, 3))
            expect_equal(dim(as.data.frame(REGION$new_data)), c(10010, 2))
        }
    )
    
    test_that(
        "KRIGING_AND_LINEAR_MODELS", {
            print("KRIGING_AND_LINEAR_MODELS:")
            REGION = LOAD_AND_DEFINE_REGION()
            MODELS = KRIGING_AND_LINEAR_MODELS(REGION)
            expect_equal(dim(MODELS$df_kriged), c(10010, 3))
            expect_equal(MODELS$model, "Matern covariance")
            expect_equal(dim(MODELS$coef_summary), c(4, 5))
        }
    )
    
    test_that(
        "PLOT_DISTRIBUTION_MAP", {
            print("PLOT_DISTRIBUTION_MAP:")
            REGION = LOAD_AND_DEFINE_REGION()
            MODELS = KRIGING_AND_LINEAR_MODELS(REGION)
            fname_map_svg = PLOT_DISTRIBUTION_MAP(REGION, MODELS)
            expect_equal(fname_map_svg, "../res/Glyphosate_resistance_distribution_SE_Australia.svg")
        }
    )
    
    test_that(
        "SAVE_SHAPEFILE", {
            print("SAVE_SHAPEFILE:")
            REGION = LOAD_AND_DEFINE_REGION()
            MODELS = KRIGING_AND_LINEAR_MODELS(REGION)
            fname_shapefile = SAVE_SHAPEFILE(MODELS)
            expect_equal(fname_shapefile, "../res/Glyphosate_resistance_gradient.shp")
        }
    )
}

# options(digits.secs=7)
# start_time = Sys.time()
# kriging_and_maps()
# end_time = Sys.time()
# print(end_time - start_time)