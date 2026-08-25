import ABC3.Found.GaloisRep.TateClassPt

/-!
# Galois (G6) 第 284 ブロック —— **★★★★★★★★★写像 `Φ : Kˣ/q^ℤ → E_q(K)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点——写像そのものができた

    Φ c := if 1 − a(c) = 0 then 0 else tatePtPair a(c) w(c) …

★★★`a(c)`・`w(c)` は類 `c` の**正規化した対**(第 213 の `exists_pair_of_class` を
`choose` したもの)。第 283 で点が類だけで決まることを示したので、
この `choose` の任意性は値に響かない。

## ★★★★★道具立てを構造体にまとめた

`TateSetup R I K` は付値 `v`、母数 `Q`(と `R` 側の `q`)、`R → K` の単射性、
`v ≥ 0` の元が `R` から来ること——を 1 つに束ねる。
★仮定が 9 本あるので、これを毎回書くと定理の頭が読めなくなる。

## ★★★場合分けは単位類でだけ切れる

`1 − a(c) = 0` になるのは `a(c) = 1`、すなわち `normRep c = 1`、すなわち **`c = 1`** の
ときだけである(`tateAOf_ne_one`)。★★★したがって

| 類 | `Φ` の値 |
|---|---|
| `c = 1` | `0`(原点) |
| `c ≠ 1` | `tatePtPair …`(**アフィン点、けっして原点でない**) |

★★これで `Φ` の単射性の「原点に行くのは単位類だけ」の段が済んだ。

## ★★次に要るもの

`Φ` の単射性(領域ごとの単射性 5 通りの貼り合わせ、第 280)、全射性、準同型性。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `TateSetup` | ★★★★★道具立ての束 |
| `tateAOf`・`tateWOf` | ★★類の正規化した対 |
| `tatePhi` | ★★★★★★★★★**写像 `Φ`** |
| `tatePhi_one` | ★★★単位類は原点に行く |
| `tatePhi_ne_zero` | ★★★★単位類でなければ原点に行かない |
| `tatePhi_eq` | ★★★★★★`Φ c` の値 |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★★★★★道具立ての束 -/

/-- ★★★★★**Tate 一意化の道具立て**——`q` と付値と `R → K` の性質をひとまとめにする。 -/
structure TateSetup (R : Type) [CommRing R] (I : Ideal R) (K : Type) [Field K]
    [Algebra R K] where
  /-- 正規化離散付値。 -/
  v : Kˣ →* Multiplicative ℤ
  /-- Tate 母数(`Kˣ` の元として)。 -/
  Q : Kˣ
  /-- `v(Q) > 0`。 -/
  hQ : 0 < vAdd v Q
  /-- Tate 母数(`R` の元として)。 -/
  q : R
  /-- `q ∈ I`。 -/
  hq : q ∈ I
  /-- 両者は一致する。 -/
  hQq : algebraMap R K q = (Q : K)
  /-- `R → K` は単射。 -/
  hinj : Function.Injective (algebraMap R K)
  /-- `v > 0` の元は `I` から来る。 -/
  hmem : ∀ x : Kˣ, 0 < vAdd v x → ∃ y ∈ I, algebraMap R K y = (x : K)
  /-- `v ≥ 0` の元は `R` から来る。 -/
  hmem0 : ∀ x : Kˣ, 0 ≤ vAdd v x → ∃ y : R, algebraMap R K y = (x : K)

/-! ## ★★類の正規化した対 -/

section Pair

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★★類 `c` の正規化した対の第 1 成分。 -/
noncomputable def tateAOf (S : TateSetup R I K) (c : Kˣ ⧸ Subgroup.zpowers S.Q) : R :=
  (exists_pair_of_class (R := R) (I := I) S.v S.Q S.hQ S.hmem S.hmem0 c).choose

/-- ★★類 `c` の正規化した対の第 2 成分。 -/
noncomputable def tateWOf (S : TateSetup R I K) (c : Kˣ ⧸ Subgroup.zpowers S.Q) : R :=
  ((exists_pair_of_class (R := R) (I := I) S.v S.Q S.hQ S.hmem S.hmem0 c).choose_spec).choose

theorem tateAOf_spec (S : TateSetup R I K) (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    algebraMap R K (tateAOf S c) = (normRep S.v S.Q S.hQ c : K) ∧ tateWOf S c ∈ I
      ∧ algebraMap R K (tateAOf S c * tateWOf S c) = (S.Q : K) :=
  ((exists_pair_of_class (R := R) (I := I) S.v S.Q S.hQ S.hmem S.hmem0 c).choose_spec).choose_spec

theorem tateWOf_mem (S : TateSetup R I K) (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    tateWOf S c ∈ I := (tateAOf_spec S c).2.1

theorem tateAOf_mul_tateWOf (S : TateSetup R I K) (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    tateAOf S c * tateWOf S c = S.q := by
  refine S.hinj ?_
  rw [(tateAOf_spec S c).2.2, S.hQq]

theorem tateAOf_one (S : TateSetup R I K) : tateAOf S 1 = 1 := by
  refine S.hinj ?_
  rw [(tateAOf_spec S 1).1, normRep_one, map_one]
  rfl

/-- ★★★★**単位類でなければ `1 − a ≠ 0`**。 -/
theorem tateAOf_ne_one (S : TateSetup R I K) {c : Kˣ ⧸ Subgroup.zpowers S.Q} (hc : c ≠ 1) :
    algebraMap R K (1 - tateAOf S c) ≠ 0 := by
  intro h
  have h0 : (1 : R) - tateAOf S c = 0 := S.hinj (by rw [h, map_zero])
  have h1 : tateAOf S c = 1 := (sub_eq_zero.1 h0).symm
  apply hc
  have h2 : (normRep S.v S.Q S.hQ c : K) = 1 := by
    rw [← (tateAOf_spec S c).1, h1, map_one]
  have h3 : normRep S.v S.Q S.hQ c = 1 := Units.ext h2
  rw [← normRep_mk S.v S.Q S.hQ c, h3]
  simp

end Pair

/-! ## ★★★★★★★★★写像 `Φ` -/

section Phi

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

/-- ★★★★★★★★★**Tate 一意化の写像** `Φ : Kˣ/q^ℤ → E_q(K)`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def tatePhi (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point :=
  if h : algebraMap R K (1 - tateAOf S c) = 0 then 0
  else tatePtPair (tateAOf S c) (tateWOf S c) S.q S.hq (tateAOf_mul_tateWOf S c)
    (isUnit_one_sub (tateWOf_mem S c)) h hΔ

/-- ★★★**単位類は原点に行く**。 -/
theorem tatePhi_one (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    tatePhi S hΔ 1 = 0 := by
  rw [tatePhi, dif_pos]
  rw [tateAOf_one, sub_self, map_zero]

/-- ★★★★★★`Φ c` の値(単位類でない場合)。 -/
theorem tatePhi_eq (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {c : Kˣ ⧸ Subgroup.zpowers S.Q} (hc : c ≠ 1) :
    tatePhi S hΔ c = tatePtPair (tateAOf S c) (tateWOf S c) S.q S.hq
      (tateAOf_mul_tateWOf S c) (isUnit_one_sub (tateWOf_mem S c))
      (tateAOf_ne_one S hc) hΔ := by
  rw [tatePhi, dif_neg (tateAOf_ne_one S hc)]

/-- ★★★★**単位類でなければ原点に行かない**。 -/
theorem tatePhi_ne_zero (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {c : Kˣ ⧸ Subgroup.zpowers S.Q} (hc : c ≠ 1) : tatePhi S hΔ c ≠ 0 := by
  rw [tatePhi, dif_neg (tateAOf_ne_one S hc), tatePtPair]
  simp

/-- ★★★★★**原点に行くのは単位類だけ**。 -/
theorem tatePhi_eq_zero_iff (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) : tatePhi S hΔ c = 0 ↔ c = 1 := by
  refine ⟨fun h => ?_, fun h => by rw [h, tatePhi_one]⟩
  by_contra hc
  exact tatePhi_ne_zero S hΔ hc h

end Phi

/-! ## ★出典の紐付け(`.src`) -/

def TateSetup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——道具立ての束)",
    sectionId := "genell-def-3-3" }

def tatePhi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——写像 Phi)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
