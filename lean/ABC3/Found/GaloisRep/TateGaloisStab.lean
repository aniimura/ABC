/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateUniformizationEquivariant

/-!
# Galois (G6) —— **★★★★★★★★★★★Tate 一意化の同変性を `Point.map` の形で**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★到達点

2026-08-27 に閉じた `tatePhi_map` は同変性を**座標の式**で述べていた
（`Φ(σ_* c) = Point.some (σ_K X) (σ_K Y)`）。
★本ファイルはそれを**群の言葉**に翻訳する:

    Point.map σ (Φ c) = Φ (σ_* c)          （`tatePhi_pointMap`）

ここで `Point.map` は **mathlib の加法群準同型**
`WeierstrassCurve.Affine.Point.map : (W⁄F).Point →+ (W⁄K).Point` である。

★★これで「`G_K` は Tate 曲線の点に作用し、`Φ` はその作用と可換」が
**posit ではなく定理**になった。

## ★★★★★★★★`σ` が `R` を固定することは自動である

`Point.map` は **`R`-代数射** `σ_A : L →ₐ[R] L` を要求する。
★`σ ∈ Gal(L/K)` は `K` を固定し、`R → K → L` なので `σ` は `R` の像を固定する
——すなわち `(σ : L →ₐ[K] L).restrictScalars R` がそのまま使える。
★★★したがって `tatePhi_map` の `σ : R →+* R` は **`RingHom.id R`** でよく、
`hcompat` は `σ_A.commutes` である。

★★★★同じ理由で `σ_U S.Q = S.Q`（母数を固定すること）も**仮説ではなく定理**になる
（`tateSetup_Q_fixed`）——`S.Q` は `algebraMap R L S.q` だからである。

## ★★★残っている仮説は 1 つだけ

`hσv : ∀ u, v(σ u) = v(u)`（付値が Galois で不変）。
★これは**付値の拡大の一意性**であり、本ファイルの守備範囲ではない。

## ★★★★★★★安定な直線の翻訳

`tate_stab_of_pointStab` は

    「点 `P = Φ(x)` の張る直線が `G_K`-安定」 ⟹ 「`σ(x) = x^c · Q^n`」

を与える。★これが `Lemma32StableLine.lean` の `lemma_3_2_i` が要求する形そのものである。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup

/-! ## ★★★★母数を固定することは自動 -/

/-- ★★★★**`R`-代数射は Tate 母数 `Q` を固定する**。

★`S.Q` の台は `algebraMap R L S.q` であり、`R`-代数射はそれを動かさない。
★★したがって `tatePhi_map` の仮説 `hσq` は**仮説ではない**。 -/
theorem tateSetup_Q_fixed {R : Type} [CommRing R] {I : Ideal R}
    {K : Type} [Field K] [Algebra R K] (S : TateSetup R I K)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K)) :
    σU S.Q = S.Q := by
  apply Units.ext
  rw [hσU, ← S.hQq, σA.commutes]

/-- ★★★**単数群への持ち上げは単射**——体の準同型が単射だからである。 -/
theorem tateSetup_σU_injective {R : Type} [CommRing R] {I : Ideal R}
    {K : Type} [Field K] [Algebra R K] (_S : TateSetup R I K)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K)) :
    Function.Injective σU := by
  intro u v h
  apply Units.ext
  have : σA (u : K) = σA (v : K) := by rw [← hσU, ← hσU, h]
  exact (σA : K →+* K).injective this

/-- ★★★**単位類でない類は `σ` で単位類に落ちない**。

★`tatePhi_map` は `σ_* c ≠ 1` を要求するので、それを `c ≠ 1` から出す段が要る。 -/
theorem tateSetup_quotMap_ne_one {R : Type} [CommRing R] {I : Ideal R}
    {K : Type} [Field K] [Algebra R K] (S : TateSetup R I K)
    (σU : Kˣ →* Kˣ) (hσq : σU S.Q = S.Q) (hinj : Function.Injective σU)
    {c : Kˣ ⧸ Subgroup.zpowers S.Q} (hc : c ≠ 1) :
    QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq) c ≠ 1 := by
  intro h
  refine hc ?_
  obtain ⟨x, rfl⟩ := QuotientGroup.mk_surjective c
  have hx : (QuotientGroup.mk (σU x) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := h
  obtain ⟨n, hn⟩ := (QuotientGroup.eq_one_iff _).1 hx
  have hn' : σU (S.Q ^ n) = σU x := by
    rw [map_zpow, hσq]; exact hn
  have hxn : x = S.Q ^ n := (hinj hn').symm
  rw [hxn]
  exact (QuotientGroup.eq_one_iff _).2 ⟨n, rfl⟩

/-! ## ★★★★★★★★★★★同変性を `Point.map` の形で -/

/-- ★★★★★★★★★★★**Tate 一意化 `Φ` は `Point.map` と可換**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

    `Point.map σ_A (Φ c) = Φ (σ_* c)`

★★左辺の `Point.map` は **mathlib の加法群準同型**である
——すなわち本定理は「`Φ` は `G_K`-加群の同型」を意味する。

★★★機構は 2 段だけ:
* 単位類は原点に行き、`Point.map` は原点を原点に送る（`map_zero`）
* それ以外では `tatePhi_map`（座標の式）に `Point.map_some` を当てる

★`rw [Point.map_some]` は `instances` 透明度で落ちる
（`(tateCurveAt q hq).map (algebraMap R K)).toAffine` と `(tateCurveAt q hq)⁄K` の食い違い）。
★★`refine (Point.map_some (S := R) σA _).trans ?_` なら通る
（`tools\lean-idioms.md`)。 -/
theorem tatePhi_pointMap {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [DecidableEq K] [Algebra R K] (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K))
    (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    Point.map (S := R) σA (tatePhi S hΔ c)
      = tatePhi S hΔ (QuotientGroup.map _ _ σU
          (zpowers_le_comap_self S.Q σU (tateSetup_Q_fixed S σA σU hσU)) c) := by
  by_cases hc : c = 1
  · subst hc
    rw [map_one, tatePhi_one]
    exact map_zero _
  · have hc' := tateSetup_quotMap_ne_one S σU (tateSetup_Q_fixed S σA σU hσU)
      (tateSetup_σU_injective S σA σU hσU) hc
    rw [tatePhi_eq S hΔ hc, tatePtPair]
    refine (Point.map_some (S := R) σA _).trans ?_
    exact (tatePhi_map S hΔ (RingHom.id R) (fun x hx => hx) (σA : K →+* K)
      (fun r => σA.commutes r) σU hσU (tateSetup_Q_fixed S σA σU hσU) hσv hc' _).symm

/-! ## ★★★★★★★安定な直線の翻訳 -/

/-- ★★★★★★★**点の側の安定性を単数の側へ移す**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

`P = Φ(x)` の張る直線が `G_K`-安定、すなわち各 `σ` について
`Point.map σ P = c·P` なら、単数の側では

    `σ(x) = x^c · Q^n`

が成り立つ。★これが `lemma_3_2_i` が要求する形そのものである。

★★入力の `Φ` は**構成したもの**でよい——`tatePhiAddEquivAll` と `fun _ => rfl` を渡せる。
★★★`Φ` を抽象のまま受けているのは、倍化の仮説（`hloc`・`hI` 等）を
本ファイルに持ち込まないためだけである。 -/
theorem tate_stab_of_pointStab {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {K L : Type} [Field K] [Field L] [DecidableEq L]
    [Algebra R K] [Algebra K L] [Algebra R L] [IsScalarTower R K L]
    (S : TateSetup R I L)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine ≠ 0)
    (Φ : Additive (Lˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R L)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (hσv : ∀ (σ : L ≃ₐ[K] L) (u : Lˣ), vAdd S.v (Units.map (σ : L →* L) u) = vAdd S.v u)
    (x : Lˣ)
    (hstab : ∀ σ : L ≃ₐ[K] L, ∃ c : ℤ,
      Point.map (S := R) ((σ : L →ₐ[K] L).restrictScalars R)
          (tatePhi S hΔ (QuotientGroup.mk x))
        = c • tatePhi S hΔ (QuotientGroup.mk x)) :
    ∀ σ : L ≃ₐ[K] L, ∃ c n : ℤ, Units.map (σ : L →* L) x = x ^ c * S.Q ^ n := by
  intro σ
  obtain ⟨c, hc⟩ := hstab σ
  have h1 : ∀ u : Lˣ, ((Units.map (σ : L →* L) u : Lˣ) : L)
      = ((σ : L →ₐ[K] L).restrictScalars R) (u : L) := fun u => rfl
  have hmap := tatePhi_pointMap S hΔ ((σ : L →ₐ[K] L).restrictScalars R)
    (Units.map (σ : L →* L)) h1 (hσv σ) (QuotientGroup.mk x)
  have key : tatePhi S hΔ (QuotientGroup.mk (Units.map (σ : L →* L) x))
      = c • tatePhi S hΔ (QuotientGroup.mk x) := by
    rw [← hc, hmap]; rfl
  have e1 : Φ (Additive.ofMul ((QuotientGroup.mk x : Lˣ ⧸ Subgroup.zpowers S.Q) ^ c))
      = c • Φ (Additive.ofMul (QuotientGroup.mk x)) := by
    rw [← map_zsmul]
    rfl
  have h2 : Φ (Additive.ofMul (QuotientGroup.mk (Units.map (σ : L →* L) x)))
      = Φ (Additive.ofMul ((QuotientGroup.mk x : Lˣ ⧸ Subgroup.zpowers S.Q) ^ c)) := by
    rw [hΦ, key, e1, hΦ]
  have h3 : (QuotientGroup.mk (Units.map (σ : L →* L) x) : Lˣ ⧸ Subgroup.zpowers S.Q)
      = QuotientGroup.mk (x ^ c) := Additive.ofMul.injective (Φ.injective h2)
  have h4 : (Units.map (σ : L →* L) x) * (x ^ c)⁻¹ ∈ Subgroup.zpowers S.Q := by
    rw [← QuotientGroup.eq_one_iff, QuotientGroup.mk_mul, h3, QuotientGroup.mk_inv,
      mul_inv_cancel]
  obtain ⟨n, hn⟩ := h4
  refine ⟨c, n, ?_⟩
  have hn' : S.Q ^ n = (Units.map (σ : L →* L) x) * (x ^ c)⁻¹ := hn
  rw [hn', ← mul_assoc, mul_comm (x ^ c) _, mul_assoc, mul_inv_cancel, mul_one]

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`(Tate 一意化)の同変性を、消費側(`Lemma 3.2, (i)`)が使う形へ
翻訳する段である。 -/

def tateSetup_Q_fixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(R-代数射は Tate 母数を固定すること)",
    sectionId := "genell-def-3-3" }

def tatePhi_pointMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化は Point.map と可換——同変性の群の形)",
    sectionId := "genell-def-3-3" }

def tate_stab_of_pointStab.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(点の側の安定性を単数の側へ移す段)",
    sectionId := "genell-lemma-3-2" }

def tatePhi_pointMap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tatePhi_map(同変性——座標の式、2026-08-27)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhi_map") 15,
    .citation "[mathlib]" "WeierstrassCurve.Affine.Point.map(点の加法群準同型)"
      (.inMathlib "WeierstrassCurve.Affine.Point.map") 15,
    .implicitStep
      ("★仮説として残るのは hσv(付値が Galois で不変)だけである" ++
       "——これは付値の拡大の一意性であり、本ファイルの守備範囲ではない") 15 ]

end ABC3.Found.GaloisRep
