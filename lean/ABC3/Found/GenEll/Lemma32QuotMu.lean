/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Lemma32StableLine

/-!
# [GenEll] Lemma 3.2, (ii) —— **`E/μ_l` の Tate 母数は `q_E^l`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★一意化の側では 1 行の群論になる

`E = K̄ˣ/q^ℤ`（Tate 一意化、`Found/GaloisRep/TateCurveWitness.lean`）とすると:

* `μ_l ⊂ K̄ˣ` の像が原文の **`μ_l ⊂ E`** である（`mu_inj` が埋め込みであることを言う）
* `E/μ_l = K̄ˣ/(μ_l·q^ℤ)`
* `K̄ˣ` は代数閉なので `x ↦ x^l` は**全射で核が `μ_l`**
* その像で `q^ℤ` は `(q^l)^ℤ` に移る

したがって

    E/μ_l ≅ K̄ˣ/(q^l)^ℤ

——★★**すなわち `q_{E′} = q_E^l`** である（`quotMuEquiv`）。

★★★そして `v_K(q^l) = l·v_K(q)`（`vAdd_pow`）なので、原文の

`deg_∞(E′) = l·deg_∞(E)`

が出る。★`deg_∞ = v_K(q_E)·log #(𝓞_K/𝔪)` は `Definition 3.3` で取ってある。

## ★★★★核になった等式

    ker(Kˣ --x↦x^l--> Kˣ --mk--> Kˣ/(q^l)^ℤ) = μ_l ⊔ q^ℤ

これが `ker_pow_mk_eq` である。★`⊆` は `x^l = (q^l)^m ⟹ (x·q^{-m})^l = 1`、
`⊇` は両方の生成元を送るだけ。★★あとは準同型定理と第 3 同型定理を継ぐ。

## ★★逸脱の記録（CLAUDE.md の「逸脱」）

★**原文は `μ_l ⊂ E` を「有限平坦部分群スキーム」として、`E′ = E/μ_l` を
スキームの商として作る。** 本ファイルは**一意化の側（`K̄ˣ/q^ℤ`）の群の商**で語る。
★★両者を繋ぐのは Tate 一意化（`uniformization_of_split`、建っている）だが、
その同型を商と両立させる段は本ファイルには入っていない。

★★★**`x ↦ x^l` の全射性を仮定として受ける**（原文の `K̄` は代数閉なので成り立つ）。
`Subgroup.zpowers` の側は `q` が 1 の冪根でないことだけを使う。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep Subgroup QuotientGroup

/-! ## ★★★★核の等式 -/

/-- ★★★★**`ker(x ↦ x^l mod (q^l)^ℤ) = μ_l ⊔ q^ℤ`**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★`⊆`: `x^l = (q^l)^m` なら `(x·q^{-m})^l = 1` なので `x ∈ μ_l·q^ℤ`。
★★`⊇`: `μ_l` の元は `l` 乗が `1`、`q` は `l` 乗が `q^l`。どちらも核に入る。 -/
theorem ker_pow_mk_eq {K : Type} [Field K] (l : ℕ) (q : Kˣ) :
    ((QuotientGroup.mk' (Subgroup.zpowers (q ^ l))).comp (powMonoidHom l)).ker
      = (powMonoidHom l : Kˣ →* Kˣ).ker ⊔ Subgroup.zpowers q := by
  have happ : ∀ y : Kˣ, (powMonoidHom l : Kˣ →* Kˣ) y = y ^ l := fun _ => rfl
  apply le_antisymm
  · intro x hx
    simp only [MonoidHom.mem_ker, MonoidHom.coe_comp, Function.comp_apply,
      QuotientGroup.mk'_apply, QuotientGroup.eq_one_iff, happ] at hx
    obtain ⟨m, hm⟩ := Subgroup.mem_zpowers_iff.1 hx
    have h1 : (x : Kˣ) ^ l = q ^ ((l:ℤ) * m) := by
      rw [← hm, ← zpow_natCast q l, ← zpow_mul]
    have hxq : (x * (q ^ m)⁻¹) ^ l = 1 := by
      rw [mul_pow, h1, ← zpow_natCast ((q ^ m)⁻¹) l, inv_zpow, ← zpow_mul, ← zpow_neg,
        ← zpow_add]
      simp [mul_comm]
    refine Subgroup.mem_sup.2 ⟨x * (q ^ m)⁻¹, ?_, q ^ m, Subgroup.zpow_mem_zpowers q m, by group⟩
    rw [MonoidHom.mem_ker, happ]
    exact hxq
  · apply sup_le
    · intro x hx
      rw [MonoidHom.mem_ker, happ] at hx
      simp only [MonoidHom.mem_ker, MonoidHom.coe_comp, Function.comp_apply,
        QuotientGroup.mk'_apply, QuotientGroup.eq_one_iff, happ, hx]
      exact one_mem _
    · intro x hx
      obtain ⟨m, hm⟩ := Subgroup.mem_zpowers_iff.1 hx
      simp only [MonoidHom.mem_ker, MonoidHom.coe_comp, Function.comp_apply,
        QuotientGroup.mk'_apply, QuotientGroup.eq_one_iff, happ]
      rw [← hm, ← zpow_natCast (q ^ m) l, ← zpow_mul]
      refine Subgroup.mem_zpowers_iff.2 ⟨m, ?_⟩
      rw [← zpow_natCast q l, ← zpow_mul, mul_comm (l:ℤ) m]

/-! ## ★★★★★★★★`Lemma 3.2, (ii)` —— `q_{E′} = q_E^l` -/

/-- ★★★★★★★★**[GenEll] Lemma 3.2, (ii)** —— `E/μ_l` は母数 `q^l` の Tate 曲線。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★これが原文の `q_{E′} = q_E^l` である。★機構は 3 段:

1. `ker_pow_mk_eq` —— `ker(x ↦ x^l mod (q^l)^ℤ) = μ_l ⊔ q^ℤ`
2. `quotientQuotientEquivQuotient`（第 3 同型定理）—— `(Kˣ/q^ℤ)/(μ_l の像) ≅ Kˣ/(μ_l ⊔ q^ℤ)`
3. `quotientKerEquivOfSurjective`（準同型定理）—— `Kˣ/ker ≅ Kˣ/(q^l)^ℤ`

★★`hsurj`（`x ↦ x^l` が全射）は原文の `K̄` が代数閉であることに当たる。 -/
noncomputable def quotMuEquiv {K : Type} [Field K] (l : ℕ) (q : Kˣ)
    (hsurj : Function.Surjective (powMonoidHom l : Kˣ →* Kˣ)) :
    ((Kˣ ⧸ Subgroup.zpowers q) ⧸
        ((powMonoidHom l : Kˣ →* Kˣ).ker.map (QuotientGroup.mk' (Subgroup.zpowers q))))
      ≃* (Kˣ ⧸ Subgroup.zpowers (q ^ l)) := by
  set H : Subgroup Kˣ := (powMonoidHom l : Kˣ →* Kˣ).ker with hH
  set N : Subgroup Kˣ := H ⊔ Subgroup.zpowers q with hN
  have hmapsup : N.map (QuotientGroup.mk' (Subgroup.zpowers q))
      = H.map (QuotientGroup.mk' (Subgroup.zpowers q)) := by
    rw [hN, Subgroup.map_sup]
    simp only [sup_eq_left]
    intro y hy
    obtain ⟨z, hz, rfl⟩ := Subgroup.mem_map.1 hy
    obtain ⟨m, rfl⟩ := Subgroup.mem_zpowers_iff.1 hz
    refine Subgroup.mem_map.2 ⟨1, one_mem _, ?_⟩
    simp only [QuotientGroup.mk'_apply, map_one]
    exact ((QuotientGroup.eq_one_iff _).2 (Subgroup.zpow_mem_zpowers q m)).symm
  refine (QuotientGroup.quotientMulEquivOfEq (by rw [← hmapsup])).trans ?_
  refine (QuotientGroup.quotientQuotientEquivQuotient (Subgroup.zpowers q) N le_sup_right).trans ?_
  refine (QuotientGroup.quotientMulEquivOfEq (by rw [hN, ← ker_pow_mk_eq l q])).trans ?_
  exact QuotientGroup.quotientKerEquivOfSurjective _ (by
    intro y
    obtain ⟨z, hz⟩ := QuotientGroup.mk_surjective (s := Subgroup.zpowers (q ^ l)) y
    obtain ⟨x, hx⟩ := hsurj z
    exact ⟨x, by simpa [hx] using hz⟩)

/-- ★★★★**`v_K(q^l) = l·v_K(q)`** —— 原文の `deg_∞(E′) = l·deg_∞(E)`。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★`deg_∞(E) = v_K(q_E)·log #(𝓞_K/𝔪)`（`Definition 3.3`）なので、
本補題と `quotMuEquiv` を継ぐと原文の等式になる。 -/
theorem vAdd_pow {K : Type} [Field K] (v : Kˣ →* Multiplicative ℤ) (l : ℕ) (q : Kˣ) :
    vAdd v (q ^ l) = l * vAdd v q := by
  show Multiplicative.toAdd (v (q ^ l)) = _
  rw [map_pow]
  simp [vAdd]

/-! ## ★★★★★★G2 非空虚／G3 負の対照 -/

/-- ★★★★**`μ_l` は `E` に埋め込まれる** —— `q` が 1 の冪根でないことを使う。

★★これが無いと「`μ_l ⊂ E` を割る」という原文の操作が空虚になりうる
（`μ_l` の像が `1` に潰れていたら `E/μ_l = E` で `q_{E′} = q_E`）。
★★★**`q` が 1 の冪根でないことが効いている**のがここである。 -/
theorem mu_inj {K : Type} [Field K] {l : ℕ} (hl : 0 < l) (q : Kˣ)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0) :
    ∀ x ∈ (powMonoidHom l : Kˣ →* Kˣ).ker,
      (QuotientGroup.mk' (Subgroup.zpowers q)) x = 1 → x = 1 := by
  intro x hx hx1
  have hxl : x ^ l = 1 := hx
  rw [QuotientGroup.mk'_apply, QuotientGroup.eq_one_iff] at hx1
  obtain ⟨m, hm⟩ := Subgroup.mem_zpowers_iff.1 hx1
  have hq0 : q ^ ((l:ℤ) * m) = 1 := by
    rw [mul_comm, zpow_mul, hm, zpow_natCast, hxl]
  have hm0 : (l:ℤ) * m = 0 := hqinf _ hq0
  have hmz : m = 0 := by
    have hlne : (l:ℤ) ≠ 0 := by exact_mod_cast hl.ne'
    exact (mul_eq_zero.1 hm0).resolve_left hlne
  rw [← hm, hmz, zpow_zero]

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** 本ファイルと `Lemma32StableLine.lean` を
合わせても、**一意化を `G_K`-同変に、また商と両立するように持ち上げる段**が
残っている（原文はスキームの水準で語る）。 -/

def ker_pow_mk_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(核の等式 ker(x ↦ x^l mod (q^l)^ℤ) = μ_l ⊔ q^ℤ)",
    sectionId := "genell-lemma-3-2" }

def quotMuEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(E/μ_l は母数 q^l の Tate 曲線——一意化の側)",
    sectionId := "genell-lemma-3-2" }

def quotMuEquiv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ker_pow_mk_eq(核の等式)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_pow_mk_eq") 15,
    .citation "[mathlib]" "QuotientGroup.quotientQuotientEquivQuotient(第 3 同型定理)"
      (.inMathlib "QuotientGroup.quotientQuotientEquivQuotient") 15,
    .citation "[mathlib]" "QuotientGroup.quotientKerEquivOfSurjective(準同型定理)"
      (.inMathlib "QuotientGroup.quotientKerEquivOfSurjective") 15,
    .implicitStep
      ("★逸脱: 原文は μ_l ⊂ E を有限平坦部分群スキームとして、E′ = E/μ_l を" ++
       "スキームの商として作る。本ファイルは一意化の側(K̄ˣ/q^ℤ)の群の商で語る。" ++
       "両者を繋ぐ Tate 一意化は建っているが、商と両立させる段は入っていない") 15,
    .implicitStep
      ("★★hsurj(x ↦ x^l が全射)は原文の K̄ が代数閉であることに当たる。" ++
       "Subgroup.zpowers の側は q が 1 の冪根でないことだけを使う") 15 ]

def vAdd_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(v_K(q^l) = l·v_K(q)——deg_∞(E′) = l·deg_∞(E) の実質)",
    sectionId := "genell-lemma-3-2" }

def mu_inj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(★μ_l が E に埋め込まれること——空虚封じ)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GenEll
