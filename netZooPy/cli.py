#!/usr/bin/python3

import click

from netZooPy import command_line as cl

@click.group()
def cli():
    pass

cli.add_command(cl.panda)
cli.add_command(cl.lioness)
cli.add_command(cl.condor)
cli.add_command(cl.bonobo)
cli.add_command(cl.otterlioness)
cli.add_command(cl.otter)