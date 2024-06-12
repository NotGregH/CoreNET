#!/usr/bin/env python

"""
This script sets up the database
for the VPlant 2.0 networking tool.

Version 1.2.0:
    - Ready to run in command line
    - Moved indexing to data insert
    - Got rid of the reverse interaction table
"""


__author__ = 'Gregory Hamilton'
__version__ = '1.2.0'
__license__ = 'MIT'
__email__ = "gah324@nyu.edu"

############
# Modules
############

import argparse as argp
from sqlalchemy import Column, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine


############
# DB setup - ORM
############

Base = declarative_base()


class Nodes(Base):

    __tablename__ = 'Nodes'
    node_id = Column(Integer, primary_key=True, autoincrement=True)
    node_name = Column(String(100), unique=True, nullable=False)


class Edges(Base):

    __tablename__ = 'Edges'
    edge_id = Column(Integer, primary_key=True, autoincrement=True)
    edge_name = Column(String(100), nullable=False, unique=True)


class Interactions(Base):

    __tablename__ = 'Interactions'
    id = Column(Integer, primary_key=True, autoincrement=True)
    node_1_id = Column(Integer, ForeignKey('Nodes.node_id'), nullable=False)
    node_2_id = Column(Integer, ForeignKey('Nodes.node_id'), nullable=False)
    edge_id = Column(Integer, ForeignKey('Edges.edge_id'), nullable=False)


############
# Functions
############


def create_db(db_name):

    engine = create_engine('sqlite:///'+db_name)
    Base.metadata.create_all(engine)

    return


############
############

if __name__ == '__main__':
    parser = argp.ArgumentParser()
    parser.add_argument('-d','--dbname',
                        help='The name you want to assign to the database file')
    args = parser.parse_args()
    create_db(args.dbname)
