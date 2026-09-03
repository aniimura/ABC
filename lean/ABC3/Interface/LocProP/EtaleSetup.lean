import Mathlib.Algebra.Group.Hom.Defs

/-!
# [LocProP] §0 —— étale cohomology の posit(`Interface`)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、物理 p.14-15。

`Definition 0.3`(arithmetic first Chern class)・`Lemma 0.4`(étale cohomology と
群 cohomology の比較)は、スキーム上の**層の étale cohomology**(`H^i(X_K, F)`)と
**Kummer 完全列**の連結射を要求する。★★**mathlib にはスキームの étale site も
étale cohomology も 0 件**(2026-09-04、`grep -i "etale.*cohomology" -r Mathlib`
相当を `.cache/mathlib-index.txt` で実測、`groupCohomology`(有限群/profinite群の
群コホモロジー)しか無い)。

★これは `Found/PGC` の p 進対数や `Found/Arakelov` の直線束テンソル積と同じ形の穴
——**古典的だが mathlib に無い数学**(建設可能、だが大きい)。ゆえに `axiom` では受けず、
`structure`(データ + 述語)として posit する。 -/

namespace ABC3.Interface.LocProP

universe u v

/-- ★posit —— `Definition 0.3` が要る最小限の骨組み。

`X_K` 自体(スキーム)は posit しない——**値が要る対象**(直線束のなす群と、
そこから誘導されるコホモロジー類)だけを posit する。

★`Type*`(auto-bound universe)は `structure ... where` の中で使うと
`picAddCommGroup` 以降が「未知の識別子」になって構文解析が壊れる
(2026-09-04 実測、`tools/lean-idioms.md` に追記予定)。`universe u v` を
明示して回避する。 -/
structure EtaleSetup where
  /-- 直線束の同型類のなす群(`H¹(X_K, 𝔾_m)` に当たる)。 -/
  Pic : Type u
  picAddCommGroup : AddCommGroup Pic
  /-- `H²(X_K, ℤ_p(1))`。 -/
  H2Zp1 : Type v
  h2AddCommGroup : AddCommGroup H2Zp1
  /-- **Kummer 完全列の連結射**から誘導される `c₁ : Pic(X_K) → H²(X_K, ℤ_p(1))`。

  原文 (LocProP p.15):
  > The connecting morphism induced on ´etale cohomology by the Kummer
  > sequence then gives us a morphism H1(XK, Gm) →H2(XK, (Z/pnZ)(1)).

  ★原文は各 `n` ごとの `H²(X_K, (ℤ/p^nℤ)(1))` への射の**両立系**として `c₁` を作るが、
  本 posit はその極限だけを直接 posit する(逆極限の構成自体は要らない)。 -/
  c1 : Pic →+ H2Zp1

end ABC3.Interface.LocProP
