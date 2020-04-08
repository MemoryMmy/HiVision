mapboxgl.accessToken = 'pk.eyJ1IjoibWVtb3J5bW15IiwiYSI6ImNqODN3amlvZjAyNXkzMnIzeGYxbXR1eG8ifQ.kEPl05qH7oxxJSVp1-tAMQ';
var map = new mapboxgl.Map({
    container: 'map',
    style:  'style.json',
    center: [106, 37],
    zoom: 4
});

var cp = $("#LayerCP").colorpicker({
    format: 'rgb',
    color: 'rgba(196,39,39,0.8)'
});


map.on('load',function() {
    var layerEl = document.getElementById('layers');
    var numEl = document.getElementById('num');
    var Enter = document.getElementById('enter');

    map.addSource('test', {
        'type': 'raster',
        'tiles': ["http://www.higis.org.cn:10082/HiVision/pchinapoint/196/39/39/0.8/{z}/{x}/{y}.png"],
        //'tileSize':256,
        'minzoom': 0,
        'maxzoom': 24
    });

    map.addLayer({
        'id': 'initial',
        'type': 'raster',
        'source': 'test'
    });

    Enter.addEventListener('click', function () {
        var layerColor = cp.data('color');
        var layer_value = layerEl.value;
        var id = String(new Date().getTime());
        console.log(id);
        var layer_id = map.getStyle().layers[map.getStyle().layers.length-1].id;
        map.removeLayer(layer_id);
        map.addSource(id, {   
            'type': 'raster',
            'tiles': [`http://www.higis.org.cn:10082/HiVision/${layer_value}/${layerColor._r}/${layerColor._g}/${layerColor._b}/${layerColor._a}/{z}/{x}/{y}.png`],
          //  'tileSize':256,
            'minzoom': 0,
            'maxzoom': 24
        });
        map.addLayer({
            'id': String(new Date().getTime()),
            'type': 'raster',
            'source': id
        });
        console.log(map.getStyle().layers);
    });

    map.addControl(new mapboxgl.NavigationControl(), 'top-right');

    window.getCenterPoint = function () {
        return {
            center: map.getCenter(),
            zoom: map.getZoom() + 1   
        };
    }

    if (window.mapView) {
        map.setCenter(window.mapView.center);
        map.setZoom(window.mapView.zoom)
    }


    window.setPlace = function (data) {
        if (map.getLayer("places")) {
            map.removeLayer("places");
        }
        if (map.getSource("places")) {
            map.removeSource("places");
        }
        var features = [];
        $.each(data, function () {
            features.push({
                "type": "Feature",
                "properties": {
                    "address": this.address,
                    "lat": this.lat,
                    "lon": this.lon,
                    "name": this.address,
                    "ogc_fid": this.ogc_fid
                },
                "geometry": {
                    "type": "Point",
                    "coordinates": [this.lon, this.lat]
                }
            });
        });
        map.addLayer({
            id: "places",
            type: "circle",
            source: {
                type: "geojson",
                data: {
                    type: "FeatureCollection",
                    features: features
                }
            }
        });
        map.on('click', 'places', function (e) {
            var coordinates = e.features[0].geometry.coordinates.slice();
            var content = '<div class="higis-attr-popup" style="min-width: 400px;height: 220px;overflow: auto">';
            var feature = e.features[0];
            if (feature) {
                var props = feature.properties;
                for (var key in props) {
                    var t = key.toString().toUpperCase();
                    if (t === "IMAGE_URI") {
                        content = content + '<h4>' + t + '</h4><div><img style="width:200px; height:200px" src=' + props[key] + '></div>';
                    }
                    else if (t === "VIDEO_URI") {
                        content = content + '<h4>' + t + '</h4><video height="200" width="300" loop="loop"><source src=' + props[key] + ' type="video/mp4">您的浏览器不支持video标签</video>';
                    }
                    else {
                        content = content + '<h4>' + t + '</h4><p>' + props[key] + '</p>';
                    }
                }
            }
            content += '</div>';
            if (popup) {
                popup.remove();
            }
            while (Math.abs(e.lngLat.lng - coordinates[0]) > 180) {
                coordinates[0] += e.lngLat.lng > coordinates[0] ? 360 : -360;
            }
            popup = new mapboxgl.Popup({
                maxHeight: 200,
                minWidth: 220
            }).setLngLat(coordinates).setHTML(content).addTo(map);
        });

        // Change the cursor to a pointer when the mouse is over the places layer.
        map.on('mouseenter', 'places', function () {
            map.getCanvas().style.cursor = 'pointer';
        });

        // Change it back to a pointer when it leaves.
        map.on('mouseleave', 'places', function () {
            map.getCanvas().style.cursor = '';
        });

    }

    var popup;
    window.showAttributes = function (e) {
        var content = '<div class="higis-attr-popup" style="min-width: 400px;height: 220px;overflow: auto">';
        var feature = e.target.feature;
        var raster = e.target.raster;
        if (feature) {
            var props = feature.properties;
            for (var key in props) {
                // content += '<h4>' + key.toString().toUpperCase() + '</h4><p>' + props[key] + '</p>';
                var t = key.toString().toUpperCase();
                if (t === "IMAGE_URI") {
                    content = content + '<h4>' + t + '</h4><div><img style="width:200px; height:200px" src=' + props[key] + '></div>';
                }
                else if (t === "VIDEO_URI") {
                    content = content + '<h4>' + t + '</h4><video height="200" width="300" loop="loop"><source src=' + props[key] + ' type="video/mp4">您的浏览器不支持video标签</video>';
                }
                else {
                    content = content + '<h4>' + t + '</h4><p>' + props[key] + '</p>';
                }
            }
        } else if (raster) {
            content += '<h4>value</h4><p>' + raster.value + '</p>';
        }
        content += '</div>';
        if (popup) {
            popup.remove();
        }
        popup = new mapboxgl.Popup({
            maxHeight: 200,
            minWidth: 220
        }).setLngLat(e.latlng).setHTML(content).addTo(map);
    }
})
