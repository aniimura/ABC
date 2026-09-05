---
name: mathlib-profinite-inventory-2026-09-05
description: Ẑ は群として mathlib にある。CompactSpace Gal(Ω/F) も instance で存在。索引は無名 instance を載せない
metadata:
  type: reference
---

pin mathlib `db127794` の実測(2026-09-05)。**「無い」という過去の記録が 2 件とも誤りだった。**

## `Ẑ` は**群として**ある

- `ProfiniteGrp.ProfiniteCompletion.completion : GrpCat.{u} → ProfiniteGrp.{u}`
  (`Mathlib.Topology.Algebra.Category.ProfiniteGrp.Completion`、2026 年に追加)
- `FiniteIndexNormalSubgroup ℤ = {nℤ}` なので
  `completion (GrpCat.of (Multiplicative ℤ))` が **そのまま `Ẑ = lim ℤ/n`**。
- 付属 API: `eta` / `denseRange` / `lift` + `lift_unique`(普遍性)/ `adjunction`。
- ★`@[to_additive]` が付いていないので `Multiplicative ℤ` にする必要がある。
- **環としては無い**(`ZHat` / `ProfiniteRing` / `Ring.completion` すべて 0 件)。
  環が要るなら `AdicCompletion`(`RingTheory/AdicCompletion/`)の型テンプレを
  `Π n, ZMod n` + `ZMod.castHom` の両立でなぞって自前で書く。

→ `Found/PGC/AbsGalSurjections.lean` 冒頭の「`Ẑ` は mathlib に不在」は**群としては誤り**。

## ★`CompactSpace Gal(Ω/F)` は instance として**ある**

`Mathlib/FieldTheory/Galois/Profinite.lean:329` の**無名 instance**
(`[IsGalois k K] : CompactSpace Gal(K/k)`)。

2026-09-05 に `infer_instance` が失敗したのは**不在だからではなく**、
`IsGalois K.carrier K.closure` が文脈にインスタンスとして無かったから。
`haveI := isGalois_carrier_closure K` を先に置けば通る。
`Found/PGC/LubinTateReciprocityCompactness.lean::compactSpace_algEquiv` は
この instance を手で書き直したものなので、置換できる。

## ★★索引の穴(在庫調査の手順に効く)

`.cache/mathlib-index.txt` は**ソース解析なので無名 instance と
`@[to_additive]` 生成名を載せない**。名前で grep しても出てこない。
今回の `CompactSpace` はソース直読で初めて見つかった。
→ `Unknown constant` / `infer_instance` 失敗のとき、[[mathlib-cohomology-inventory-2026-09-05]]
と同じく「無い」と結論する前に**ソースを直読する**か、索引を無名 instance 込みに拡張すること。

## `Γ_K ↠ Ẑ` は `Ẑ` を作らずに済む(最短路)

`exists_surjective_absGal_to_unramifiedClosureGal`(`Γ_K ↠ Gal(K^ur/K)`)と
`isGalois_unramifiedClosure` が既にあるので、上の instance で
`Gal(K^ur/K)` は自動的にコンパクト・副有限。**残るのは `Gal(K^ur/K) ≅ Ẑ` 一本**。
★`exists_surjective_absGal_to_zmod` は各 `n` で独立に `∃ φ` を出すだけで
**`n` をまたいだ両立が無い**——`Ẑ` へ 1 本の射を作るならそこが真の難所。
`unramifiedClosure` 経由ならその両立は既に入っている。

## 逆極限・コンパクト性の在庫

- `IsCompact.nonempty_iInter_of_directed_nonempty_isCompact_isClosed`(**有向版**)
  ——木の ℕ 特化(`reciprocityMapLimit_surjective`)はこれに差し替えるだけで
  `(ℕ, ∣)` のような非全順序の添字に昇格できる。
- `nonempty_sections_of_finite_cofiltered_system` / `TopCat.nonempty_limitCone_of_compact_t2_cofiltered_system`
- `CategoryTheory.Limits.Types.surjective_π_app_zero_of_surjective_map`(ℕᵒᵖ 塔、有限性不要)
- `ProfiniteGrp.limitConePtAux : Subgroup (Π j, F.obj j)` は
  ABC3 の `CompatibleUnits` と**同じ作り**。
- `ABC3.Found.PGC.compatible_of_succ` は既に完全に一般(隣接両立 ⇒ 全両立)。再利用可。
- `PadicInt` は Cauchy 完備化であって逆極限ではないので、**作り方は流用できない**
  (定理側 `lift`/`lift_unique`/`ext_of_toZModPow` は手本になる)。
