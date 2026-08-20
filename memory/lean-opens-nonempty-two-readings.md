---
name: lean-opens-nonempty-two-readings
description: mathlib の `U : X.Opens` に対する `Nonempty U` は `↥U` と `↥↑U` の 2 通りに読めてしまい、instance が噛み合わない
metadata:
  type: reference
---

`U : X.Opens`（`AlgebraicGeometry`）に対して `Nonempty U` と書くと、
コアーションの経路が 2 つある：

* `↥U`（`SetLike` 経由の `{x // x ∈ U}`）—— **mathlib の
  `Scheme.germToFunctionField (U : X.Opens) [h : Nonempty U]` が要求するのはこちら**
* `↥↑U`（`Opens → Set → Sort`）—— `Set.Nonempty.to_subtype` が返すのはこちら

そのため `haveI := hU.to_subtype`（`hU : (↑U).Nonempty`）や
`haveI : Nonempty ↥U := hU.to_subtype` では
`Algebra Γ(X, U) X.functionField` の instance 探索が失敗する。

**How to apply:** `haveI : Nonempty U := ⟨⟨hU.some, hU.some_mem⟩⟩` と、
**型注釈を `Nonempty U` と素直に書いて** Lean に mathlib と同じ経路を選ばせる。
関連: [[mathlib-local-instances]]
