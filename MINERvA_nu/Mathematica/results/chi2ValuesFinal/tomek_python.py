data = np.array([[0.36184742, 0.28947794, 0.23158235, 0.18526588, 0.1482127 ,
        0.11857016, 0.1       , 0.1       , 0.1       , 0.1       ,
        0.1       , 0.1       , 0.11430328, 0.1428791 , 0.17859888,
        0.2232486 , 0.27906075, 0.34882594, 0.43603242, 0.54504053,
        0.52345074, 0.46047015, 0.55433232, 0.45230927],
       [0.        , 0.23158235, 0.18526588, 0.1482127 , 0.11857016,
        0.1       , 0.1       , 0.1       , 0.1       , 0.1       ,
        0.1       , 0.1       , 0.1       , 0.11430328, 0.1428791 ,
        0.17859888, 0.2232486 , 0.27906075, 0.34882594, 0.43603242,
        0.46047015, 0.36837612, 0.46047015, 0.47850945],
       [0.        , 0.        , 0.1482127 , 0.11857016, 0.1       ,
        0.1       , 0.1       , 0.11015008, 0.11136773, 0.1       ,
        0.1       , 0.1       , 0.10193594, 0.1       , 0.11430328,
        0.1428791 , 0.17859888, 0.2232486 , 0.27906075, 0.34882594,
        0.43603242, 0.44038662, 0.47850945, 0.50891506],
       [0.        , 0.        , 0.        , 0.1       , 0.1       ,
        0.1       , 0.11015008, 0.1376876 , 0.13920966, 0.11274581,
        0.1       , 0.10193594, 0.12741993, 0.1132732 , 0.1415915 ,
        0.11430328, 0.1428791 , 0.17859888, 0.2232486 , 0.27906075,
        0.34882594, 0.43603242, 0.54504053, 0.63614383],
       [0.        , 0.        , 0.        , 0.        , 0.1       ,
        0.11015008, 0.1376876 , 0.17210949, 0.17401208, 0.14093226,
        0.11274581, 0.12741993, 0.15927491, 0.1415915 , 0.17698937,
        0.1415915 , 0.16676129, 0.20845161, 0.26056451, 0.32570564,
        0.40713205, 0.50891506, 0.63614383, 0.79517978],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.1376876 , 0.17210949, 0.21513687, 0.2175151 , 0.17616533,
        0.14093226, 0.15927491, 0.19909364, 0.17698937, 0.22123671,
        0.17698937, 0.20845161, 0.26056451, 0.32570564, 0.40713205,
        0.50891506, 0.63614383, 0.79517978, 0.99397473],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.17401208, 0.2175151 , 0.27189387, 0.22020666,
        0.17616533, 0.19909364, 0.24886705, 0.22123671, 0.27654589,
        0.22123671, 0.17698937, 0.20845161, 0.26056451, 0.32570564,
        0.40713205, 0.50891506, 0.63614383, 0.79517978],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.17616533, 0.22020666, 0.27525833,
        0.22020666, 0.24886705, 0.31108381, 0.27654589, 0.34568237,
        0.27654589, 0.22123671, 0.25431671, 0.31789589, 0.39736986,
        0.49671233, 0.62089041, 0.77611302, 0.97014127],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.17616533, 0.22020666,
        0.27525833, 0.23239165, 0.29048956, 0.34568237, 0.43210296,
        0.34568237, 0.27654589, 0.31789589, 0.39736986, 0.49671233,
        0.62089041, 0.77611302, 0.97014127, 1.21267659],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.17698937,
        0.22123671, 0.27654589, 0.34568237, 0.43210296, 0.5401287 ,
        0.43210296, 0.34568237, 0.39736986, 0.49671233, 0.62089041,
        0.77611302, 0.97014127, 1.21267659, 1.51584573],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.27654589, 0.34568237, 0.43210296, 0.5401287 , 0.67516087,
        0.5401287 , 0.43210296, 0.34568237, 0.39736986, 0.49671233,
        0.62089041, 0.77611302, 0.97014127, 1.21267659],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.27654589, 0.34568237, 0.43210296, 0.5401287 ,
        0.54861309, 0.43889047, 0.35111238, 0.40267959, 0.50334948,
        0.62918685, 0.78648357, 0.98310446, 1.22888057],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.4282785 , 0.53534812, 0.43210296,
        0.46125773, 0.51383926, 0.41107141, 0.32885713, 0.40267959,
        0.50334948, 0.62918685, 0.78648357, 0.98310446],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.4282785 , 0.53534812,
        0.50944035, 0.63680044, 0.51383926, 0.41107141, 0.32885713,
        0.40267959, 0.50334948, 0.62918685, 0.78648357],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.50944035,
        0.63680044, 0.51383926, 0.64229907, 0.51383926, 0.41107141,
        0.32885713, 0.40267959, 0.50334948, 0.62918685],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.79600055, 0.64229907, 0.51383926, 0.41107141, 0.32885713,
        0.2630857 , 0.32214367, 0.40267959, 0.50334948],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.51383926, 0.41107141, 0.32885713, 0.29785995,
        0.23828796, 0.25771493, 0.32214367, 0.40267959],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.44761355, 0.35809084, 0.28647267,
        0.29619982, 0.23695986, 0.25771493, 0.32214367],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.44761355, 0.35809084,
        0.28647267, 0.29619982, 0.23695986, 0.25771493],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.28647267,
        0.22917814, 0.23695986, 0.29619982, 0.23695986],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.28647267, 0.22917814, 0.23695986, 0.29619982],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.28647267, 0.25339101, 0.27906075],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.31673876, 0.34882594],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.43603242]])
