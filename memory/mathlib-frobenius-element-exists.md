---
name: mathlib-frobenius-element-exists
description: Frobenius 元は mathlib にある(IsArithFrobAt / arithFrobAt)。名前が違うので旧 regex に当たらなかった
metadata:
  type: reference
---

**`Mathlib.RingTheory.Frobenius`(Andrew Yang, 2025)に Frobenius 元がある。**

- `IsArithFrobAt (R) {S} [CommRing R] [CommRing S] [Algebra R S] {M} [Monoid M]
  [MulSemiringAction M S] [SMulCommClass M R S] : M → Ideal S → Prop`
- `IsArithFrobAt.exists_of_isInvariant` / `arithFrobAt (R) (G) (Q) : G`
- `isConj_arithFrobAt : Ideal.under R Q = Ideal.under R Q' →
  IsConj (arithFrobAt R G Q) (arithFrobAt R G Q')`
- `AlgHom.IsArithFrobAt` + 15 本の補題(`eq_of_isUnramifiedAt`(一意性)ほか)

★**流儀が `Gal(L/K)` ではない**。`G` が `S` に作用し `R` が固定環という
`RingTheory/Invariant` の形。要るインスタンスは `[Algebra.IsInvariant R S G]`・
`[SMulCommClass G R S]`・`[Finite (S ⧸ Q)]`。共役類は `IsConj` で受ける。

**Why:** `Skeleton/NumberField/Chebotarev.lean` の `.needs` は
「`FrobeniusElement` で grep、0 件(2026-08-20)」と記録していたが**誤り**だった。
名前が `IsArithFrobAt` なので旧 regex に当たらない。
[[mathlib-cohomology-inventory-2026-09-05]] と同じ「不在の誤り」で、
2026-09-05〜06 の 2 日で **5 件目**(`ULift.field` / `continuousCohomology` / `Ẑ` /
`CompactSpace Gal` / これ)。

**How to apply:** 「無い」と書く前に `node tools/absent-recheck.mjs --try '<正規表現>'`
で件数を測る(第 1019 で `check.mjs` の G11 が新規に要求するようになった)。
★**名前で引けないときは概念の別名を複数試す**。今回は `Frobenius` 単独なら当たった。

なお依然 0 件(2026-09-06 再測): `Chebotarev|Tchebotarev` / `DirichletDensity` /
`ArtinReciprocity|artinMap` / `RayClassGroup` / `SplitsCompletely` /
`decompositionGroup` / `analyticDensity|naturalDensity|primeDensity`。
