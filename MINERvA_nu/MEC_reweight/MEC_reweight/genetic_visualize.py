import numpy as np
from matplotlib import pyplot


@np.vectorize
def function(x, y):
    x = x / 10 - 3
    y = y / 10 - 3
    return np.cos(x**2 + y**2) / (np.abs(x**2 + y**2) + 1)


def new_population(old):
    points = []
    scores = []

    for x, y in old:
        points.append(x)
        scores.append(y)

    reduction = min(scores)

    scores = [s - reduction for s in scores]


def mutate(point):
    new_point = point.copy()
    new_point += np.random.normal(size=2) * 10
    new_point = np.clip(new_point, a_min=0, a_max=100)
    return new_point


population = [np.random.random((2,)) * 100 for _ in range(10)]

population_with_score = [(p, function(*p)) for p in population]
population_with_score.sort(key=lambda x: x[1], reverse=True)
sorted_population = np.array([x[0] * 10 for x in population_with_score])

x = y = np.linspace(0, 100, 1000)
xv, yv = np.meshgrid(x, y)
z = function(xv, yv)

pyplot.imshow(z)
pyplot.savefig('function.pdf')
pyplot.show()


pyplot.imshow(z)
pyplot.scatter(sorted_population.T[0], sorted_population.T[1], c='green', s=5)

for _ in range(100):
    new_population = [population_with_score[0][0]]
    for _ in range(len(population) - 2):
        new_population.append(mutate(new_population[0]))
    new_population.append(np.random.random((2,)) * 100)

    population = new_population

    population_with_score = [(p, function(*p)) for p in population]
    population_with_score.sort(key=lambda x: x[1], reverse=True)
    sorted_population = np.array([x[0] * 10 for x in population_with_score])
    pyplot.scatter(sorted_population.T[0],
                   sorted_population.T[1], c='black', s=0.1)

pyplot.scatter(sorted_population.T[0], sorted_population.T[1], c='red', s=10)
pyplot.savefig('dots.pdf')
pyplot.show()
