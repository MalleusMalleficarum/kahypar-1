/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/
#pragma once

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "kahypar/io/sql_plottools_serializer.h"
#include "kahypar/partition/evolutionary/edge_frequency.h"
#include "kahypar/partition/evolutionary/population.h"


namespace kahypar {
namespace combine {
static constexpr bool debug = false;

Individual partitions(Hypergraph& hg,
                      const Parents& parents,
                      Context& context) {
  DBG << V(context.evolutionary.action.decision());
  DBG << "Parent 1: initial" << V(parents.first.fitness());
  DBG << "Parent 2: initial" << V(parents.second.fitness());
  context.evolutionary.parent1 = &parents.first.partition();
  context.evolutionary.parent2 = &parents.second.partition();
#ifndef NDEBUG
  hg.setPartition(parents.first.partition());
  ASSERT(parents.first.fitness() == metrics::km1(hg));
  DBG << "initial" << V(metrics::km1(hg)) << V(metrics::imbalance(hg, context));
  hg.reset();
  if (!context.evolutionary.action.requires().invalidation_of_second_partition) {
    hg.setPartition(parents.second.partition());
    ASSERT(parents.second.fitness() == metrics::km1(hg));
    DBG << "initial" << V(metrics::km1(hg)) << V(metrics::imbalance(hg, context));
    hg.reset();
  }
#endif
  hg.reset();
  HypernodeID original_contraction_limit_multiplier = context.coarsening.contraction_limit_multiplier;
  if (context.evolutionary.unlimited_coarsening_contraction) {
    context.coarsening.contraction_limit_multiplier = 1;
  }
  Partitioner().partition(hg, context);
  context.coarsening.contraction_limit_multiplier = original_contraction_limit_multiplier;
  DBG << "Offspring" << V(metrics::km1(hg)) << V(metrics::imbalance(hg, context));
  ASSERT(metrics::km1(hg) <= std::min(parents.first.fitness(), parents.second.fitness()));
  io::serializer::serializeEvolutionary(context, hg);
  io::printEvolutionaryInformation(context);
  Individual indi = Individual(hg);
  return Individual(hg);
}


Individual usingTournamentSelection(Hypergraph& hg, const Context& context, const Population& population) {
  Context temporary_context(context);

  temporary_context.evolutionary.action =
    Action { meta::Int2Type<static_cast<int>(EvoDecision::combine)>() };
  temporary_context.coarsening.rating.rating_function = RatingFunction::heavy_edge;
  temporary_context.coarsening.rating.partition_policy = RatingPartitionPolicy::evolutionary;

  const auto& parents = population.tournamentSelect();


  return combine::partitions(hg, parents, temporary_context);
}



Individual edgeFrequency(Hypergraph& hg, const Context& context, const Population& population) {
  hg.reset();
  Context temporary_context(context);

  temporary_context.evolutionary.action =
    Action { meta::Int2Type<static_cast<int>(EvoDecision::combine)>(),
             meta::Int2Type<static_cast<int>(EvoCombineStrategy::edge_frequency)>() };

  temporary_context.coarsening.rating.rating_function = RatingFunction::edge_frequency;
  temporary_context.coarsening.rating.partition_policy = RatingPartitionPolicy::normal;
  temporary_context.coarsening.rating.heavy_node_penalty_policy =
    HeavyNodePenaltyPolicy::edge_frequency_penalty;

  temporary_context.evolutionary.edge_frequency =
    computeEdgeFrequency(population.listOfBest(context.evolutionary.edge_frequency_amount),
                         hg.initialNumEdges());

  DBG << V(temporary_context.evolutionary.action.decision());

  Partitioner().partition(hg, temporary_context);
  DBG << "final result" << V(metrics::km1(hg)) << V(metrics::imbalance(hg, context));
  io::serializer::serializeEvolutionary(temporary_context, hg);
  io::printEvolutionaryInformation(temporary_context);
  return Individual(hg);
}

Individual usingTournamentSelectionAndEdgeFrequency(Hypergraph& hg,
                                                    const Context& context,
                                                    const Population& population) {
  Context temporary_context(context);

  temporary_context.evolutionary.action =
    Action { meta::Int2Type<static_cast<int>(EvoDecision::combine)>() };

  temporary_context.coarsening.rating.rating_function = RatingFunction::edge_frequency;
  temporary_context.coarsening.rating.partition_policy = RatingPartitionPolicy::evolutionary;
  temporary_context.evolutionary.edge_frequency =
    computeEdgeFrequency(population.listOfBest(context.evolutionary.edge_frequency_amount),
                         hg.initialNumEdges());

  const auto& parents = population.tournamentSelect();

  // TODO remove

  // temporary_context.evolutionary.parent1 = &parents.first.get().partition();
  // temporary_context.evolutionary.parent2 = &parents.second.get().partition();

  DBG << V(temporary_context.evolutionary.action.decision());
  return combine::partitions(hg, parents, temporary_context);
}


// TODO(andre) is this even viable?

}  // namespace combine
}  // namespace kahypar
