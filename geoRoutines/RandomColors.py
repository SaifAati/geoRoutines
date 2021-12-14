# Author : Saif Aati
# Contact: SAIF AATI  <saif@caltech.edu> <saifaati@gmail.com>
# Copyright (C) 2020
#https://gist.github.com/adewes/5884820

import random

def get_random_color(pastel_factor=0.5):
    return [(x + pastel_factor) / (1.0 + pastel_factor) for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]]


def color_distance(c1, c2):
    return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])


def generate_new_color(existing_colors, pastel_factor=0.5):
    max_distance = None
    best_color = None
    for i in range(0, 100):
        color = get_random_color(pastel_factor=pastel_factor)
        if not existing_colors:
            return color
        best_distance = min([color_distance(color, c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color

def GenerateColors (N,pastel_factor=0.9):
    colors = []

    for i in range(0, N):
        colors.append(generate_new_color(colors, pastel_factor=0.9))

    # print("Your colors:", colors)
    return colors

if __name__ == '__main__':

    # To make your color choice reproducible, uncomment the following line:
    # random.seed(10)
    GenerateColors(N=10)