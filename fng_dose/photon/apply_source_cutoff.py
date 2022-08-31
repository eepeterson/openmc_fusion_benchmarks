import openmc
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('cutoff', type=float)
args = parser.parse_args()

settings = openmc.Settings.from_xml()
for src in settings.source:
    dist = src.energy
    mask = dist.x > args.cutoff
    dist.x = dist.x[mask]
    dist.p = dist.p[mask]
    src.strength = dist.p.sum()
settings.export_to_xml()
