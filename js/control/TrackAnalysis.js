/**
 * Provides track analysis functionality.
 *
 * Takes the detailed way tags from brouter-server's response
 * and creates tables with distributions of way types, surfaces,
 * and smoothness values.
 *
 * On hovering/click a table row the corresponding track segments
 * are highlighted on the map.
 *
 * @type {L.Class}
 */
BR.TrackAnalysis = L.Class.extend({
    /**
     * @type {Object}
     */
    options: {
        overlayStyle: {
            color: 'yellow',
            opacity: 0.8,
            weight: 8,
            // show above quality coding (pane defined in RoutingPathQuality.js)
            pane: 'routingQualityPane',
        },
    },

    /**
     * The total distance of the whole track, recalculate on each `update()` call.
     *
     * @type {float}
     */
    totalRouteDistance: 0.0,

    /**
     * @param {Map} map
     * @param {object} options
     */
    initialize(map, options) {
        this.map = map;
        L.setOptions(this, options);
    },

    /**
     * @type {?BR.TrackEdges}
     */
    trackEdges: null,

    /**
     * @type {?L.Polyline}
     */
    trackPolyline: null,

    /**
     * true when tab is shown, false when hidden
     *
     * @type {boolean}
     */
    active: false,

    /**
     * Called by BR.Sidebar when tab is activated
     */
    show() {
        this.active = true;
        this.options.requestUpdate(this);
    },

    /**
     * Called by BR.Sidebar when tab is deactivated
     */
    hide() {
        this.active = false;
    },

    /**
     * Everytime the track changes this method is called:
     *
     * - calculate statistics (way type, max speed, surface, smoothness)
     *   for the whole track
     * - renders statistics tables
     * - create event listeners which allow to hover/click a
     *   table row for highlighting matching track segments
     *
     * @param {Polyline} polyline
     * @param {Array} segments - route segments between waypoints
     */
    update(polyline, segments) {
        if (!this.active) {
            return;
        }

        if (segments.length === 0) {
            $('#track_statistics').html('');
            if (this.highlightedSegments) {
                this.map.removeLayer(this.highlightedSegments);
                this.highlightedSegments = null;
            }
            if (this.highlightedSegment) {
                this.map.removeLayer(this.highlightedSegment);
                this.highlightedSegment = null;
            }
            return;
        }

        this.trackPolyline = polyline;
        this.trackEdges = new BR.TrackEdges(segments);

        const analysis = this.calcStats(polyline, segments);

        this.render(analysis);

        $('.track-analysis-table tr').hover(L.bind(this.handleHover, this), L.bind(this.handleHoverOut, this));
        $('.track-analysis-table tbody').on('click', 'tr', L.bind(this.toggleSelected, this));
    },

    /**
     * This method calculates the radius of curvature between three points
     */
    calculateRadiusOfCurvature(point1, point2, point3) {
        // const latLongPoints = [
        //     [point1[1], point1[0]],
        //     [point2[1], point2[0]],
        //     [point3[1], point3[0]]
        // ];

        const measuredPointLatitude = point2[1];
        const measuredPointLongitude = point2[0];
        const originPointLatitude = point1[1];
        const originPointLongitude = point1[0];
        const destinationPointLatitude = point3[1];
        const destinationPointLongitude = point3[0];

        // console.log("[", originPointLatitude, ",", originPointLongitude, "],[", measuredPointLatitude, ",", measuredPointLongitude, "],[", destinationPointLatitude, ",", destinationPointLongitude + "]");

        // console.log([
        //     { latitude: originPointLatitude, longitude: measuredPointLongitude },
        //     { latitude: measuredPointLatitude, longitude: measuredPointLongitude },
        //     { latitude: measuredPointLatitude, longitude: originPointLongitude },
        //     { latitude: measuredPointLatitude, longitude: measuredPointLongitude },
        //     { latitude: destinationPointLatitude, longitude: measuredPointLongitude },
        //     { latitude: measuredPointLatitude, longitude: measuredPointLongitude },
        //     { latitude: measuredPointLatitude, longitude: measuredPointLongitude },
        //     { latitude: measuredPointLatitude, longitude: destinationPointLongitude }
        // ])

        // East/west or longitude distance between points 1 and 2
        const distance1 =
            (Math.sign(measuredPointLatitude - originPointLatitude) *
                geolib.getDistance(
                    { latitude: originPointLatitude, longitude: measuredPointLongitude },
                    { latitude: measuredPointLatitude, longitude: measuredPointLongitude },
                    (accuracy = 0.1)
                )) /
            1000;

        // North/south or latitude distance between points 1 and 2
        const distance2 =
            (Math.sign(measuredPointLongitude - originPointLongitude) *
                geolib.getDistance(
                    { latitude: measuredPointLatitude, longitude: originPointLongitude },
                    { latitude: measuredPointLatitude, longitude: measuredPointLongitude },
                    (accuracy = 0.1)
                )) /
            1000;

        // East/west or longitude distance between points 2 and 3
        const distance3 =
            (Math.sign(measuredPointLatitude - destinationPointLatitude) *
                geolib.getDistance(
                    { latitude: destinationPointLatitude, longitude: measuredPointLongitude },
                    { latitude: measuredPointLatitude, longitude: measuredPointLongitude },
                    (accuracy = 0.1)
                )) /
            1000;

        // North/south or latitude distance between points 2 and 3
        const distance4 =
            (Math.sign(measuredPointLongitude - destinationPointLongitude) *
                geolib.getDistance(
                    { latitude: measuredPointLatitude, longitude: measuredPointLongitude },
                    { latitude: measuredPointLatitude, longitude: destinationPointLongitude },
                    (accuracy = 0.1)
                )) /
            1000;

        let points = [
            [distance1, distance2],
            [0, 0],
            [distance3, distance4],
        ];

        if (Math.abs(distance1) < 0.001 || Math.abs(distance3) < 0.001 || Math.abs(distance1 - distance3) < 0.001) {
            // If either of the latitudes are the same as the point we are measuring, our quadratic function will not be valid
            // Swapping x/y should not change the arc length, chord, or radius calculations.
            // We treat anything closer than a meter as being equal.
            points = [
                [distance2, distance1],
                [0, 0],
                [distance4, distance3],
            ];
        }

        // console.log(points);

        // Extract x and y coordinates
        const x = points.map((point) => point[0]);
        const y = points.map((point) => point[1]);

        // Set up the system of equations matrix A
        const A = x.map((xi) => [xi ** 2, xi, 1]);

        // Solve for the coefficients a, b, c
        const coefficients = math.lusolve(A, y);

        // Extract the coefficients
        const [a, b, c] = coefficients.map((coef) => coef[0]);

        // console.log(`a: ${a}, b: ${b}, c: ${c}`);

        // Define the derivative of the quadratic polynomial
        function derivative(x) {
            return 2 * a * x + b;
        }

        // Define the integrand for the arc length
        function integrand(x) {
            const power = math.pow(derivative(x), 2);
            const add = math.add(1, power);
            // @ts-ignore
            return math.sqrt(add);
        }

        // console.log(`Integrand: ${integrand(5)}`);

        // Define the points (x1, y1) and (x3, y3)
        const x1 = points[0][0];
        const y1 = points[0][1];
        const x3 = points[2][0];
        const y3 = points[2][1];

        // a = lower bound
        // b = upper bound
        // n = number of trapezoids
        function trapezoidalIntegral(f, a, b, n) {
            const h = (b - a) / n;
            let sum = 0.5 * (f(a) + f(b));
            for (let i = 1; i < n; i++) {
                sum += f(a + i * h);
            }
            return sum * h;
        }

        // Calculate the arc length using numerical integration
        const arcLength = trapezoidalIntegral(integrand, x1, x3, 2000);

        // Calculate the chord length
        //@ts-ignore
        const chord = math.sqrt(math.add(math.pow(x3 - x1, 2), math.pow(y3 - y1, 2)));

        if (typeof chord !== 'number') {
            throw new Error('Chord is not a number');
        }

        // console.log(`Arc length: ${arcLength}, Chord length: ${chord}`);

        function solveForRUsingSecant(a, d, r0, r1) {
            const tolerance = 1e-8;
            const maxIterations = 1000;

            let f0 = 2 * r0 * math.sin(a / (2 * r0)) - d;
            let f1 = 2 * r1 * math.sin(a / (2 * r1)) - d;

            for (let i = 0; i < maxIterations; i++) {
                if (math.abs(f1 - f0) === 0) {
                    // TODO horrible assumption, but this means the number is going to be large enough that we don't care about it.
                    return 20;
                    //throw new Error('Division by zero in secant method');
                }

                const rNext = r1 - (f1 * (r1 - r0)) / (f1 - f0);

                if (math.abs(rNext - r1) < tolerance) {
                    return rNext;
                }

                r0 = r1;
                f0 = f1;
                r1 = rNext;
                f1 = 2 * r1 * math.sin(a / (2 * r1)) - d;
            }

            throw new Error('Failed to converge for parameters ' + a + ', ' + d + ', ' + r0 + ', ' + r1);
        }

        // derivative of 2*x*math.sin(0.076/ (2 * x))-0.075
        // const radius = solveForR(arcLength, chord, chord / 2);
        const radius = solveForRUsingSecant(arcLength, chord, chord * 2, chord / 2);

        return radius;
    },

    /**
     * This method does the heavy-lifting of statistics calculation.
     *
     * What happens here?
     *
     * - loop over all route segments
     * - for each segment loop over all contained points
     * - parse and analyze the `waytags` field between two consecutive points
     * - group the values for each examined category (highway, surface, smoothness) and sum up the distances
     *   - special handling for tracks: create an entry for each tracktype (and one if the tracktype is unknown)
     * - sort the result by distance descending
     *
     * @param polyline
     * @param segments
     * @returns {Object}
     */
    calcStats(polyline, segments) {
        const analysis = {
            highway: {},
            maxspeed: {},
            surface: {},
            smoothness: {},
            curvature: {},
            grade: {},
            trainConstruction: {},
        };

        this.totalRouteDistance = 0.0;

        const myLatLngs = [];

        for (let segmentIndex = 0; segments && segmentIndex < segments.length; segmentIndex++) {
            for (
                let coordinateIndex = 1;
                coordinateIndex < segments[segmentIndex].feature.geometry.coordinates.length - 1;
                coordinateIndex++
            ) {
                if (coordinateIndex == 0) {
                    continue;
                }

                const coordinate = segments[segmentIndex].feature.geometry.coordinates[coordinateIndex];

                myLatLngs.push(L.latLng(coordinate[1], coordinate[0], coordinate[2]));
            }
        }

        var geojsonFeatures = geoDataExchange.buildGeojsonFeatures(myLatLngs, {
            interpolate: false,
            normalize: false,
        });
        var latLngToGradeMap = {};
        var latLngToCurvatureMap = {};

        for (let featureIndex = 0; featureIndex < geojsonFeatures[0]['features'].length; featureIndex++) {
            const feature = geojsonFeatures[0]['features'][featureIndex];

            for (let geoIndex = 0; geoIndex < feature['geometry']['coordinates'].length; geoIndex++) {
                latLngToGradeMap[feature['geometry']['coordinates'][geoIndex]] = feature['properties']['attributeType'];
            }
        }

        for (let segmentIndex = 0; segments && segmentIndex < segments.length; segmentIndex++) {
            for (
                let messageIndex = 1;
                messageIndex < segments[segmentIndex].feature.properties.messages.length;
                messageIndex++
            ) {
                this.totalRouteDistance += parseFloat(
                    segments[segmentIndex].feature.properties.messages[messageIndex][3]
                );
                let wayTags = segments[segmentIndex].feature.properties.messages[messageIndex][9].split(' ');
                wayTags = this.normalizeWayTags(wayTags, 'cycling');
                for (let wayTagIndex = 0; wayTagIndex < wayTags.length; wayTagIndex++) {
                    let wayTagParts = wayTags[wayTagIndex].split('=');
                    let tagName = wayTagParts[0];
                    switch (tagName) {
                        case 'highway':
                            let highwayType = wayTagParts[1];
                            let trackType = '';
                            if (highwayType === 'track') {
                                trackType = this.getTrackType(wayTags);
                                highwayType = 'Track ' + trackType;
                            }
                            if (typeof analysis.highway[highwayType] === 'undefined') {
                                analysis.highway[highwayType] = {
                                    formatted_name: i18next.t(
                                        'sidebar.analysis.data.highway.' + highwayType,
                                        highwayType
                                    ),
                                    name: wayTagParts[1],
                                    subtype: trackType,
                                    distance: 0.0,
                                };
                            }
                            analysis.highway[highwayType].distance += parseFloat(
                                segments[segmentIndex].feature.properties.messages[messageIndex][3]
                            );
                            break;
                        case 'maxspeed':
                        case 'surface':
                        case 'smoothness':
                            if (typeof analysis[tagName][wayTagParts[1]] === 'undefined') {
                                let formattedName;

                                if (tagName.indexOf('maxspeed') === 0) {
                                    formattedName = i18next.t('sidebar.analysis.data.maxspeed', {
                                        maxspeed: wayTagParts[1],
                                    });
                                } else {
                                    formattedName = i18next.t([
                                        'sidebar.analysis.data.' + tagName + '.' + wayTagParts[1],
                                        wayTagParts[1],
                                    ]);
                                }

                                analysis[tagName][wayTagParts[1]] = {
                                    formatted_name: formattedName,
                                    name: wayTagParts[1],
                                    subtype: '',
                                    distance: 0.0,
                                };
                            }
                            analysis[tagName][wayTagParts[1]].distance += parseFloat(
                                segments[segmentIndex].feature.properties.messages[messageIndex][3]
                            );
                            break;
                    }
                }
            }

            // Loop over the actual node points to determine if any violate curvature or grade rules
            for (
                let coordinateIndex = 1;
                coordinateIndex < segments[segmentIndex].feature.geometry.coordinates.length - 1;
                coordinateIndex++
            ) {
                if (coordinateIndex == 0) {
                    continue;
                }

                const coordinate = segments[segmentIndex].feature.geometry.coordinates[coordinateIndex];

                const point1 = segments[segmentIndex].feature.geometry.coordinates[coordinateIndex - 1];
                const point2 = coordinate;
                const point3 = segments[segmentIndex].feature.geometry.coordinates[coordinateIndex + 1];

                // Based on the height graph, it looks like the grade is calculated at the point we are going TO
                // So that is how we should count distance
                const grade = latLngToGradeMap[coordinate];

                const curvature = this.calculateRadiusOfCurvature(point1, point2, point3);

                latLngToCurvatureMap[coordinate] = curvature;

                if (typeof analysis.grade['valid'] === 'undefined') {
                    analysis.grade['valid'] = {
                        formatted_name: 'valid', //i18next.t('sidebar.analysis.data.trainConstruction.curvature'),
                        name: 'valid',
                        subtype: '',
                        distance: 0.0,
                    };
                }

                if (typeof analysis.grade['invalid'] === 'undefined') {
                    analysis.grade['invalid'] = {
                        formatted_name: 'invalid', //i18next.t('sidebar.analysis.data.trainConstruction.curvature'),
                        name: 'invalid',
                        subtype: '',
                        distance: 0.0,
                    };
                }

                if (typeof analysis.trainConstruction['valid'] === 'undefined') {
                    analysis.trainConstruction['valid'] = {
                        formatted_name: 'valid', //i18next.t('sidebar.analysis.data.trainConstruction.curvature'),
                        name: 'valid',
                        subtype: '',
                        distance: 0.0,
                    };
                }

                if (typeof analysis.trainConstruction['invalid'] === 'undefined') {
                    analysis.trainConstruction['invalid'] = {
                        formatted_name: 'invalid', //i18next.t('sidebar.analysis.data.trainConstruction.curvature'),
                        name: 'invalid',
                        subtype: '',
                        distance: 0.0,
                    };
                }

                if (curvature > 0.4) {
                    analysis.trainConstruction['valid']['distance'] += geolib.getDistance(
                        { longitude: point1[0], latitude: point1[1] },
                        { longitude: point2[0], latitude: point2[1] }
                    );
                } else {
                    analysis.trainConstruction['invalid']['distance'] += geolib.getDistance(
                        { longitude: point1[0], latitude: point1[1] },
                        { longitude: point2[0], latitude: point2[1] }
                    );
                }

                if (grade > 2 || grade < -2) {
                    analysis.grade['invalid']['distance'] += geolib.getDistance(
                        { longitude: point1[0], latitude: point1[1] },
                        { longitude: point2[0], latitude: point2[1] }
                    );
                } else {
                    analysis.grade['valid']['distance'] += geolib.getDistance(
                        { longitude: point1[0], latitude: point1[1] },
                        { longitude: point2[0], latitude: point2[1] }
                    );
                }
            }

            const trackLatLngs = this.trackPolyline.getLatLngs();

            for (let i = 0; i < trackLatLngs.length; i++) {
                trackLatLngs[i]['grade'] =
                    latLngToGradeMap[`${trackLatLngs[i]['lng']},${trackLatLngs[i]['lat']},${trackLatLngs[i]['alt']}`];
                trackLatLngs[i]['curvature'] =
                    latLngToCurvatureMap[
                        `${trackLatLngs[i]['lng']},${trackLatLngs[i]['lat']},${trackLatLngs[i]['alt']}`
                    ];
            }
        }

        BR.Grade = latLngToGradeMap;
        BR.Curvature = latLngToCurvatureMap;

        return this.sortAnalysisData(analysis);
    },

    /**
     * Normalize the tag name.
     *
     * Motivation: The `surface` and `smoothness` tags come in different variations,
     * e.g. `surface`, `cycleway:surface` etc. We're only interested
     * in the tag which matches the given routing type. All other variations
     * are dropped. If no specialized surface/smoothness tag is found, the default value
     * is returned, i.e. `smoothness` or `surface`.
     *
     * Also, maxspeed comes in different variations, e.g. `maxspeed`, `maxspeed:forward`,
     * `maxspeed:backward`. Depending on the existence of the `reversedirection` field
     * we can select the correct value.
     *
     * @param wayTags - tags + values for a way segment
     * @param routingType - currently only 'cycling' is supported, can be extended in the future (walking, driving, etc.)
     * @returns {*[]}
     */
    normalizeWayTags(wayTags, routingType) {
        let normalizedWayTags = {};
        let surfaceTags = {};
        let smoothnessTags = {};
        for (let wayTagIndex = 0; wayTagIndex < wayTags.length; wayTagIndex++) {
            let wayTagParts = wayTags[wayTagIndex].split('=');
            const tagName = wayTagParts[0];
            const tagValue = wayTagParts[1];

            if (tagName === 'surface') {
                surfaceTags.default = tagValue;
                continue;
            }
            if (tagName.indexOf(':surface') !== -1) {
                let tagNameParts = tagName.split(':');
                surfaceTags[tagNameParts[0]] = tagValue;
                continue;
            }

            if (tagName === 'smoothness') {
                smoothnessTags.default = tagValue;
                continue;
            }
            if (tagName.indexOf(':smoothness') !== -1) {
                let tagNameParts = tagName.split(':');
                smoothnessTags[tagNameParts[0]] = tagValue;
                continue;
            }

            if (tagName === 'maxspeed:forward' && !wayTags.includes('reversedirection=yes')) {
                normalizedWayTags['maxspeed'] = tagValue;
                continue;
            }
            if (tagName === 'maxspeed:backward' && wayTags.includes('reversedirection=yes')) {
                normalizedWayTags['maxspeed'] = tagValue;
                continue;
            }
            if (tagName === 'maxspeed') {
                normalizedWayTags[tagName] = tagValue;
                continue;
            }

            normalizedWayTags[tagName] = tagValue;
        }

        switch (routingType) {
            case 'cycling':
                if (typeof surfaceTags.cycleway === 'string') {
                    normalizedWayTags.surface = surfaceTags.cycleway;
                } else if (typeof surfaceTags.default === 'string') {
                    normalizedWayTags.surface = surfaceTags.default;
                }
                if (typeof smoothnessTags.cycleway === 'string') {
                    normalizedWayTags.smoothness = smoothnessTags.cycleway;
                } else if (typeof smoothnessTags.default === 'string') {
                    normalizedWayTags.smoothness = smoothnessTags.default;
                }
                break;
            default:
                if (typeof surfaceTags.default === 'string') {
                    normalizedWayTags.surface = surfaceTags.default;
                }
                if (typeof smoothnessTags.default === 'string') {
                    normalizedWayTags.smoothness = smoothnessTags.default;
                }
        }

        return this.wayTagsToArray(normalizedWayTags);
    },

    /**
     * Transform analysis data for each type into an array, sort it
     * by distance descending and convert it back to an object.
     *
     * @param {Object} analysis
     *
     * @returns {Object}
     */
    sortAnalysisData(analysis) {
        const analysisSortable = {};
        const result = {};

        for (const type in analysis) {
            if (!analysis.hasOwnProperty(type)) {
                continue;
            }

            result[type] = {};
            analysisSortable[type] = [];

            for (const name in analysis[type]) {
                if (!analysis[type].hasOwnProperty(name)) {
                    continue;
                }
                analysisSortable[type].push(analysis[type][name]);
            }

            if (type === 'maxspeed') {
                analysisSortable[type].sort(function (a, b) {
                    return parseInt(a.name) - parseInt(b.name);
                });
            } else {
                analysisSortable[type].sort(function (a, b) {
                    return b.distance - a.distance;
                });
            }

            for (let j = 0; j < analysisSortable[type].length; j++) {
                result[type][analysisSortable[type][j].formatted_name] = analysisSortable[type][j];
            }
        }

        return result;
    },

    /**
     * Extract the tracktype from a waytags string.
     * If no tracktype is found 'unknown' is returned.
     *
     * @param {string[]} wayTags
     * @returns {string}
     */
    getTrackType(wayTags) {
        for (let i = 0; i < wayTags.length; i++) {
            const wayTagParts = wayTags[i].split('=');
            if (wayTagParts[0] === 'tracktype') {
                return wayTagParts[1];
            }
        }

        return 'unknown';
    },

    /**
     * @param {Object} analysis
     */
    render(analysis) {
        const $content = $('#track_statistics');

        $content.html('');
        $content.append($(`<h4 class="track-analysis-heading">${i18next.t('sidebar.analysis.header.highway')}</h4>`));
        $content.append(this.renderTable('highway', analysis.highway));
        $content.append($(`<h4 class="track-analysis-heading">${i18next.t('sidebar.analysis.header.surface')}</h4>`));
        $content.append(this.renderTable('surface', analysis.surface));
        $content.append(
            $(`<h4 class="track-analysis-heading">${i18next.t('sidebar.analysis.header.smoothness')}</h4>`)
        );
        $content.append(this.renderTable('smoothness', analysis.smoothness));
        $content.append($(`<h4 class="track-analysis-heading">${i18next.t('sidebar.analysis.header.maxspeed')}</h4>`));
        $content.append(this.renderTable('maxspeed', analysis.maxspeed));
        $content.append($(`<h4 class="track-analysis-heading">${i18next.t('sidebar.analysis.header.grade')}</h4>`));
        $content.append(this.renderTable('grade', analysis.grade));
        $content.append(
            $(`<h4 class="track-analysis-heading">${i18next.t('sidebar.analysis.header.trainconstruction')}</h4>`)
        );
        $content.append(this.renderTable('trainconstruction', analysis.trainConstruction));
    },

    /**
     * Renders an analysis table.
     *
     * @param {string} type
     * @param {Array} data
     * @returns {jQuery}
     */
    renderTable(type, data) {
        let index;
        const $table = $(`<table data-type="${type}" class="mini stripe dataTable track-analysis-table"></table>`);
        const $thead = $('<thead></thead>');
        $thead.append(
            $('<tr>')
                .append(
                    `<th class="track-analysis-header-category">${i18next.t('sidebar.analysis.table.category')}</th>`
                )
                .append(
                    $(`<th class="track-analysis-header-distance">${i18next.t('sidebar.analysis.table.length')}</th>`)
                )
        );
        $table.append($thead);
        const $tbody = $('<tbody></tbody>');

        let totalDistance = 0.0;

        for (index in data) {
            if (!data.hasOwnProperty(index)) {
                continue;
            }
            const $row = $(`<tr data-name="${data[index].name}" \
                data-subtype="${data[index].subtype}" \
                data-distance="${data[index].distance}"></tr>`);
            $row.append(`<td class="track-analysis-title">${data[index].formatted_name}</td>`);
            $row.append(`<td class="track-analysis-distance">${this.formatDistance(data[index].distance)} km</td>`);
            $tbody.append($row);
            totalDistance += data[index].distance;
        }

        if (totalDistance < this.totalRouteDistance) {
            $tbody.append(
                $(`<tr data-name="internal-unknown" data-distance="${this.totalRouteDistance - totalDistance}"></tr>`)
                    .append($(`<td class="track-analysis-title">${i18next.t('sidebar.analysis.table.unknown')}</td>`))
                    .append(
                        $(
                            `<td class="track-analysis-distance">${this.formatDistance(
                                this.totalRouteDistance - totalDistance
                            )} km</td>`
                        )
                    )
            );
        }

        $table.append($tbody);

        $table.append(
            $('<tfoot></tfoot>')
                .append('<tr></tr>')
                .append($(`<td>${i18next.t('sidebar.analysis.table.total_known')}</td>`))
                .append(
                    $(
                        `<td class="track-analysis-distance track-analysis-distance-total">${this.formatDistance(
                            totalDistance
                        )} km</td>`
                    )
                )
        );

        return $table;
    },

    /**
     * Format a distance with two decimal places.
     *
     * @param {number} meters
     * @returns {string}
     */
    formatDistance(meters) {
        return (meters / 1000).toLocaleString(undefined, { minimumFractionDigits: 2, maximumFractionDigits: 2 });
    },

    handleHover(event) {
        const $tableRow = $(event.currentTarget);
        const $table = $tableRow.parents('table').first();
        const dataType = $table.data('type');
        const dataName = $tableRow.data('name');
        const trackType = $tableRow.data('subtype');

        const polylinesForDataType = this.getPolylinesForDataType(dataType, dataName, trackType);

        this.highlightedSegments = L.layerGroup(polylinesForDataType).addTo(this.map);
    },

    handleHoverOut() {
        this.map.removeLayer(this.highlightedSegments);
    },

    toggleSelected(event) {
        const tableRow = event.currentTarget;
        const $table = $(tableRow).parents('table').first();
        const dataType = $table.data('type');
        const dataName = $(tableRow).data('name');
        const trackType = $(tableRow).data('subtype');

        if (tableRow.classList.toggle('selected')) {
            if (this.highlightedSegment) {
                this.map.removeLayer(this.highlightedSegment);
                this.selectedTableRow.classList.remove('selected');
            }
            this.highlightedSegment = L.layerGroup(this.getPolylinesForDataType(dataType, dataName, trackType)).addTo(
                this.map
            );
            this.selectedTableRow = tableRow;

            return;
        }

        this.map.removeLayer(this.highlightedSegment);
        this.selectedTableRow = null;
        this.highlightedSegment = null;
    },

    /**
     * Searching each track edge if it matches the requested
     * arguments (type, name, subtype if type == track). If the
     * track edge matches the search, create a Leaflet polyline
     * and add it to the result array.
     *
     * @param {string} dataType - `highway`, `surface`, `smoothness`
     * @param {string} dataName - `primary`, `track, `asphalt`, etc.
     * @param {string} trackType - the tracktype is passed here (e.g.
     * `grade3`), but only in the case that `dataName` is `track`
     *
     * @returns {Polyline[]}
     */
    getPolylinesForDataType(dataType, dataName, trackType) {
        const polylines = [];
        const trackLatLngs = this.trackPolyline.getLatLngs();

        if (dataType === 'trainconstruction' || dataType === 'grade') {
            const curPolyline = [];
            let isValid = false;
            for (let i = 0; i < trackLatLngs.length; i++) {
                const curvature = trackLatLngs[i].curvature;
                const grade = trackLatLngs[i].grade;

                if (grade === undefined || curvature === undefined) {
                    console.log('No grade or curvature found for this point');
                    continue;
                }

                isValid = false;

                if (dataType === 'trainconstruction' && curvature > 0.4) {
                    isValid = true;
                } else if (dataType === 'grade' && grade < 2) {
                    isValid = true;
                }

                if ((dataName === 'valid' && isValid) || (dataName === 'invalid' && !isValid)) {
                    curPolyline.push(trackLatLngs[i]);
                } else {
                    if (curPolyline.length > 1) {
                        polylines.push(L.polyline(curPolyline, this.options.overlayStyle));
                    }
                    curPolyline.length = 0;
                    curPolyline.push(trackLatLngs[i]);
                }
            }

            if (curPolyline.length > 1 && ((dataName === 'valid' && isValid) || (dataName === 'invalid' && !isValid))) {
                polylines.push(L.polyline(curPolyline, this.options.overlayStyle));
                curPolyline.length = 0;
            }
        } else {
            for (let i = 0; i < this.trackEdges.edges.length; i++) {
                if (this.wayTagsMatchesData(trackLatLngs[this.trackEdges.edges[i]], dataType, dataName, trackType)) {
                    const matchedEdgeIndexStart = i > 0 ? this.trackEdges.edges[i - 1] : 0;
                    const matchedEdgeIndexEnd = this.trackEdges.edges[i] + 1;
                    polylines.push(
                        L.polyline(
                            trackLatLngs.slice(matchedEdgeIndexStart, matchedEdgeIndexEnd),
                            this.options.overlayStyle
                        )
                    );
                }
            }
        }

        return polylines;
    },

    /**
     * Examine the way tags string if it matches the data arguments.
     * Special handling for implicit defined dataName 'internal-unknown'
     * which matches if a tag-pair is missing. Special handling for
     * tracktypes again.
     *
     * @param {string} wayTags - The way tags as provided by brouter, e.g.
     * `highway=secondary surface=asphalt smoothness=good`
     * @param {string} dataType - `highway`, `surface`, `smoothness`
     * @param {string} dataName - `primary`, `track, `asphalt`, etc.
     * @param {string} trackType - the tracktype is passed here (e.g.
     * `grade3`), but only in the case that `dataName` is `track`
     *
     * @returns {boolean}
     */
    wayTagsMatchesData(wayTags, dataType, dataName, trackType) {
        const parsed = this.wayTagsToObject(wayTags);

        switch (dataType) {
            case 'highway':
                if (dataName === 'track') {
                    if (trackType === 'unknown' && parsed.highway === 'track' && !parsed.tracktype) {
                        return true;
                    }

                    return typeof parsed.tracktype === 'string' && parsed.tracktype === trackType;
                } else if (dataName === 'internal-unknown' && typeof parsed.highway !== 'string') {
                    return true;
                }

                return typeof parsed.highway === 'string' && parsed.highway === dataName;
            case 'surface':
                return this.singleWayTagMatchesData('surface', parsed, dataName);
            case 'smoothness':
                return this.singleWayTagMatchesData('smoothness', parsed, dataName);
            case 'maxspeed':
                return this.singleWayTagMatchesData('maxspeed', parsed, dataName);
        }

        return false;
    },

    singleWayTagMatchesData(category, parsedData, lookupValue) {
        if (typeof lookupValue === 'number') {
            lookupValue = lookupValue.toString();
        }

        let foundValue = null;

        // We need to handle `maxspeed:forward` and `maxspeed:backward` separately
        // from all other tags, because we need to consider the `reversedirection`
        // tag.
        // Test URL: http://localhost:3000/#map=15/52.2292/13.6204/standard&lonlats=13.61948,52.231611;13.611327,52.227431
        if (
            category === 'maxspeed' &&
            parsedData.hasOwnProperty('maxspeed:forward') &&
            !parsedData.hasOwnProperty('reversedirection')
        ) {
            foundValue = parsedData['maxspeed:forward'];
        }
        if (
            category === 'maxspeed' &&
            parsedData.hasOwnProperty('maxspeed:backward') &&
            parsedData.hasOwnProperty('reversedirection') &&
            parsedData.reversedirection === 'yes'
        ) {
            foundValue = parsedData['maxspeed:backward'];
        }

        // if the special handling for `maxspeed` didn't find a result,
        // check wayTags for matching property:
        if (foundValue === null && parsedData.hasOwnProperty(category)) {
            foundValue = parsedData[category];
        }

        if (lookupValue === 'internal-unknown' && foundValue === null) {
            return true;
        }

        return foundValue === lookupValue;
    },

    /**
     * Transform a way tags string into an object, for example:
     *
     * 'highway=primary surface=asphalt' => { highway: 'primary', surface: 'asphalt' }
     *
     * @param wayTags - The way tags as provided by brouter, e.g.
     * `highway=secondary surface=asphalt smoothness=good`
     *
     * @returns {object}
     */
    wayTagsToObject(wayTags) {
        let result = {};
        const wayTagPairs = wayTags.feature.wayTags.split(' ');

        for (let j = 0; j < wayTagPairs.length; j++) {
            const wayTagParts = wayTagPairs[j].split('=');
            result[wayTagParts[0]] = wayTagParts[1];
        }

        return result;
    },

    /**
     * Transform a way tags object into an array representation, for example:
     *
     * { 'highway' : 'path', 'surface' : 'sand' } => ['highway=path', 'surface=sand']
     *
     * @param wayTags - The way tags in object representation
     *
     * @returns {object}
     */
    wayTagsToArray(wayTags) {
        let wayTagsArray = [];
        for (let wayTagKey in wayTags) {
            wayTagsArray.push(wayTagKey + '=' + wayTags[wayTagKey]);
        }

        return wayTagsArray;
    },

    /**
     * @param point1 - The first point
     * @param point2 - The second point
     * @param point3 - The third point
     *
     * @returns {number} - The curvature of the three points
     */
    calculateCurvature(point1, point2, point3) {
        // TODO: This is almost certainly wrong. I'm pretty sure this needs to be the arc radius approach.
        const a = Math.sqrt(Math.pow(point2[0] - point1[0], 2) + Math.pow(point2[1] - point1[1], 2));
        const b = Math.sqrt(Math.pow(point3[0] - point2[0], 2) + Math.pow(point3[1] - point2[1], 2));
        const c = Math.sqrt(Math.pow(point1[0] - point3[0], 2) + Math.pow(point1[1] - point3[1], 2));

        return a + b + c;
    },
});
