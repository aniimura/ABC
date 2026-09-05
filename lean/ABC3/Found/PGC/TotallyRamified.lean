import ABC3.Found.PGC.UnramifiedSubextension

/-!
# 完全分岐と、不分岐との交わり

`UnramifiedSubextension.lean` で `K(x) ≤ K^ur ⟺ K(x)/K が不分岐` まで来た。
相互律の全体像(`Γ_K^ab ≅ (K^×)^∧`)に必要な**線型無関係性**
「完全分岐 ∩ K^ur = K」へ向けて、まずその両端を用意する。

* `IsTotallyRamifiedAdjoin K x := inertiaDegree K x = 1`(慣性次数が 1)
* `finrank_eq_one_of_isUnramified_of_isTotallyRamified`——
  不分岐かつ完全分岐なら次数 1(`e·f = [K(x):K]` の直接の帰結)
* `isUnit_adjoinIntegers_iff`——`adjoinIntegers K x` の単数はノルム 1 の元。
  部分拡大に沿った整数環の包含が局所準同型であることを示すのに要る。
* `inertiaDegree_le_of_adjoin_le`——慣性次数の単調性(剰余体は伸びる一方)
* **`isTotallyRamified_of_le`**——完全分岐は部分拡大に遺伝する
* **`finrank_eq_one_of_mem_unramifiedClosure_of_le`**——「完全分岐 ∩ K^ur = K」

## ★配管の記録: 中間体をまたぐ `rfl` は kernel を止める

`K(x) ≤ K(y)` に沿った整数環の包含 `adjoinIntegers K x → adjoinIntegers K y` は
数学的には自明(中間体のノルムは `K.closure` のノルムの制限そのもの)だが、
**Lean では素直に書けない**:

* `IntermediateField.inclusion hle` を使うと、`((inclusion hle z : ↥K⟮y⟯) : K.closure)
  = ((z : ↥K⟮x⟯) : K.closure)` の `rfl` が **kernel deterministic timeout**(実測 60 秒)。
* `⟨(z : K.closure), hle z.2⟩` と素直に書けば、**中間体の元 1 層なら 0.06 秒**で通る
  (`‖⟨(z : K.closure), hle w.2⟩‖ = ‖w‖` は `rfl`)。
* ところが `adjoinIntegers`(中間体の部分環)の元は**2 層**なので、
  同じ形の `def` でも kernel が落ちる(実測 60 秒)。

つまり「1 層は速い・2 層は落ちる」。

**★回避策が効いた**: `def` の中で membership を defeq に頼らせず、
ノルム保存 `norm_mk_of_le` を**先に**補題として用意して `rw` で渡すと、
60 秒の kernel timeout が **0.12 秒**になる。これで包含
`adjoinIntegersIncl` / `adjoinIntegersRingHom` が作れ、剰余体の比較
`𝓀_x ↪ 𝓀_y` から**慣性次数の単調性**が出て、目標の

**`finrank_eq_one_of_mem_unramifiedClosure_of_le`
(= 完全分岐拡大の中の不分岐な部分は自明 = 「完全分岐 ∩ K^ur = K」)**

まで届いた。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- **完全分岐**——慣性次数が 1(剰余体が伸びない)。 -/
def IsTotallyRamifiedAdjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] : Prop :=
  inertiaDegree K x = 1

/-- **不分岐かつ完全分岐なら次数 1**——`e·f = [K(x):K]` の直接の帰結。
「完全分岐 ∩ K^ur = K」の**片方の端**。 -/
theorem finrank_eq_one_of_isUnramified_of_isTotallyRamified (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (hu : IsUnramifiedAdjoin K x) (ht : IsTotallyRamifiedAdjoin K x) :
    Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = 1 := by
  have h := ramificationIndex_mul_inertiaDegree K x
  rw [show ramificationIndex K x = 1 from hu, show inertiaDegree K x = 1 from ht] at h
  omega

/-- `x ∈ K^ur` かつ `K(x)/K` が完全分岐なら `K(x) = K`。 -/
theorem finrank_eq_one_of_mem_unramifiedClosure_of_isTotallyRamified (K : PAdicLocalField p)
    (x : K.closure) (hx : x ∈ unramifiedClosure K) (ht : IsTotallyRamifiedAdjoin K x) :
    Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = 1 :=
  finrank_eq_one_of_isUnramified_of_isTotallyRamified K x
    ((mem_unramifiedClosure_iff_isUnramified K x).mp hx) ht

/-- **`adjoinIntegers K x` の単数はノルム 1 の元**。
部分拡大に沿った整数環の包含が局所準同型であることを示すのに要る。 -/
theorem isUnit_adjoinIntegers_iff (K : PAdicLocalField p) (x : K.closure)
    (z : adjoinIntegers K x) :
    IsUnit z ↔ ‖(z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ = 1 := by
  constructor
  · rintro ⟨u, rfl⟩
    have h1 : ((u : adjoinIntegers K x) : IntermediateField.adjoin K.carrier
        ({x} : Set K.closure)) * ((↑u⁻¹ : adjoinIntegers K x) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = 1 :=
      congrArg Subtype.val u.mul_inv
    have h2 := congrArg norm h1
    rw [norm_mul, norm_one] at h2
    have hu : ‖((u : adjoinIntegers K x) : IntermediateField.adjoin K.carrier
      ({x} : Set K.closure))‖ ≤ 1 := (u : adjoinIntegers K x).2
    have hv : ‖((↑u⁻¹ : adjoinIntegers K x) : IntermediateField.adjoin K.carrier
      ({x} : Set K.closure))‖ ≤ 1 := (↑u⁻¹ : adjoinIntegers K x).2
    nlinarith [norm_nonneg ((u : adjoinIntegers K x) : IntermediateField.adjoin K.carrier
      ({x} : Set K.closure)), norm_nonneg ((↑u⁻¹ : adjoinIntegers K x) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  · intro hz
    have hz0 : ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ≠ 0 := by
      intro h0
      rw [h0, norm_zero] at hz
      exact zero_ne_one hz
    refine ⟨⟨z, ⟨((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))⁻¹, ?_⟩, ?_, ?_⟩,
      rfl⟩
    · show ‖((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))⁻¹‖ ≤ 1
      rw [norm_inv, hz, inv_one]
    · apply Subtype.ext
      show ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
        * ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))⁻¹ = 1
      exact mul_inv_cancel₀ hz0
    · apply Subtype.ext
      show ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))⁻¹
        * ((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))) = 1
      exact inv_mul_cancel₀ hz0

/-! ## 部分拡大に沿った整数環の包含

★上の docstring に書いた壁の**回避策が効いた**: `def` の中で membership を
defeq に頼らせず、ノルム保存 `norm_mk_of_le` を**先に**補題として用意して
`rw` で渡すと、60 秒の kernel timeout が **0.12 秒**になる。 -/

/-- 中間体をまたぐノルムの保存(1 層なので `rfl` が速い)。 -/
theorem norm_mk_of_le (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (w : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ‖(⟨(w : K.closure), hle w.2⟩ : IntermediateField.adjoin K.carrier
      ({y} : Set K.closure))‖ = ‖w‖ := rfl

/-- `K(x) ≤ K(y)` に沿った整数環の包含(写像)。 -/
noncomputable def adjoinIntegersIncl (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (z : adjoinIntegers K x) : adjoinIntegers K y :=
  ⟨⟨((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure),
      hle (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2⟩,
    by
      show ‖(⟨((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure),
        hle (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2⟩ :
          IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ ≤ 1
      rw [norm_mk_of_le K hle]
      exact z.2⟩

@[simp] theorem norm_adjoinIntegersIncl (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (z : adjoinIntegers K x) :
    ‖((adjoinIntegersIncl K hle z : adjoinIntegers K y)
      : IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖
      = ‖((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)))‖ :=
  norm_mk_of_le K hle _

/-- 環準同型としての包含。 -/
noncomputable def adjoinIntegersRingHom (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    adjoinIntegers K x →+* adjoinIntegers K y where
  toFun := adjoinIntegersIncl K hle
  map_one' := by apply Subtype.ext; apply Subtype.ext; rfl
  map_mul' _ _ := by apply Subtype.ext; apply Subtype.ext; rfl
  map_zero' := by apply Subtype.ext; apply Subtype.ext; rfl
  map_add' _ _ := by apply Subtype.ext; apply Subtype.ext; rfl

@[simp] theorem adjoinIntegersRingHom_apply (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (z : adjoinIntegers K x) :
    adjoinIntegersRingHom K hle z = adjoinIntegersIncl K hle z := rfl

/-- 包含は局所準同型(ノルムを保つから)。 -/
instance isLocalHom_adjoinIntegersRingHom (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    IsLocalHom (adjoinIntegersRingHom K hle) := by
  refine ⟨fun z hz => ?_⟩
  rw [isUnit_adjoinIntegers_iff] at hz ⊢
  rw [adjoinIntegersRingHom_apply, norm_adjoinIntegersIncl] at hz
  exact hz

/-! ## 慣性次数の単調性、そして「完全分岐 ∩ K^ur = K」 -/

/-- **剰余体は伸びる一方**——`K(x) ≤ K(y)` なら `𝓀_x ↪ 𝓀_y`。 -/
theorem card_residueField_le_of_adjoin_le (K : PAdicLocalField p) {x y : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    residueDegree K x ≤ residueDegree K y := by
  haveI := isLocalRing_adjoinIntegers K x
  haveI := isLocalRing_adjoinIntegers K y
  exact Nat.card_le_card_of_injective _
    (IsLocalRing.ResidueField.map (adjoinIntegersRingHom K hle)).injective

/-- **慣性次数は単調**。 -/
theorem inertiaDegree_le_of_adjoin_le (K : PAdicLocalField p) {x y : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    inertiaDegree K x ≤ inertiaDegree K y := by
  have h := card_residueField_le_of_adjoin_le K hle
  rw [residueDegree_eq_residueCard_pow K x, residueDegree_eq_residueCard_pow K y] at h
  have hq2 : 1 < Nat.card 𝓀[K.carrier] := by
    haveI : Fintype 𝓀[K.carrier] := Fintype.ofFinite _
    rw [Nat.card_eq_fintype_card]
    exact Fintype.one_lt_card
  exact (Nat.pow_le_pow_iff_right hq2).mp h

/-- **★★★★★完全分岐は部分拡大に遺伝する**。 -/
theorem isTotallyRamified_of_le (K : PAdicLocalField p) {x y : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (hty : IsTotallyRamifiedAdjoin K y) : IsTotallyRamifiedAdjoin K x := by
  have hle' := inertiaDegree_le_of_adjoin_le K hle
  rw [show inertiaDegree K y = 1 from hty] at hle'
  have hpos : 0 < inertiaDegree K x := by
    have h := ramificationIndex_mul_inertiaDegree K x
    have hn : 0 < Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := Module.finrank_pos
    rcases Nat.eq_zero_or_pos (inertiaDegree K x) with h0 | h0
    · rw [h0, Nat.mul_zero] at h; omega
    · exact h0
  show inertiaDegree K x = 1
  omega

/-- **★★★★★★「完全分岐 ∩ K^ur = K」**——完全分岐拡大の中の不分岐な部分は自明。
相互律の全体像(`Γ_K^ab ≅ (K^×)^∧`)に必要な線型無関係性。 -/
theorem finrank_eq_one_of_mem_unramifiedClosure_of_le (K : PAdicLocalField p) {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (hx : x ∈ unramifiedClosure K) (hty : IsTotallyRamifiedAdjoin K y) :
    Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = 1 :=
  finrank_eq_one_of_isUnramified_of_isTotallyRamified K x
    ((mem_unramifiedClosure_iff_isUnramified K x).mp hx)
    (isTotallyRamified_of_le K hle hty)

/-! ## 完全分岐の判定条件 -/

/-- **完全分岐 ⟺ 剰余体が伸びない**——`residueDegree = q`。
`isUnramifiedAdjoin_iff_residueDegree`(不分岐 ⟺ `residueDegree = q^{[K(x):K]}`)の双対。 -/
theorem isTotallyRamifiedAdjoin_iff_residueDegree (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    IsTotallyRamifiedAdjoin K x ↔ residueDegree K x = Nat.card 𝓀[K.carrier] := by
  rw [residueDegree_eq_residueCard_pow K x]
  have hq2 : 2 ≤ Nat.card 𝓀[K.carrier] := by
    haveI : Fintype 𝓀[K.carrier] := Fintype.ofFinite _
    have h1 : 1 < Fintype.card 𝓀[K.carrier] := Fintype.one_lt_card
    rw [Nat.card_eq_fintype_card]
    omega
  constructor
  · intro h; rw [show inertiaDegree K x = 1 from h, pow_one]
  · intro h
    nth_rewrite 2 [← pow_one (Nat.card 𝓀[K.carrier])] at h
    exact Nat.pow_right_injective hq2 h

/-! ## Newton 多角形の一段: Eisenstein の根のノルム

★一般の超距離ノルム体で成り立つ補題。次の段(Lubin-Tate 塔が完全分岐である
ことの証明)で使う: `‖α‖^n = ‖a_0‖ = ‖π‖` から値群の指数が `n` 以上になり、
`e ≤ e·f = n` と合わせて `e = n`・`f = 1` が出る。 -/

/-- モニック多項式の根で、低次係数がすべて `‖a_0‖ < 1` 以下なら、根は開単位球に入る。 -/
theorem norm_lt_one_of_monic_root {F : Type*} [NormedField F] [IsUltrametricDist F]
    {n : ℕ} (_hn : 0 < n) (a : ℕ → F) {α : F}
    (hroot : α ^ n + ∑ i ∈ Finset.range n, a i * α ^ i = 0)
    (hle : ∀ i, i < n → ‖a i‖ ≤ ‖a 0‖) (hlt : ‖a 0‖ < 1) (hne : a 0 ≠ 0) :
    ‖α‖ < 1 := by
  by_contra hge
  rw [not_lt] at hge
  have hα0 : α ≠ 0 := by
    intro h0
    rw [h0, norm_zero] at hge
    linarith
  have hαpos : (0:ℝ) < ‖α‖ ^ n := by positivity
  have hbound : ‖∑ i ∈ Finset.range n, a i * α ^ i‖ ≤ ‖a 0‖ * ‖α‖ ^ n := by
    refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
    intro i hi
    rw [Finset.mem_range] at hi
    rw [norm_mul, norm_pow]
    refine mul_le_mul (hle i hi) ?_ (by positivity) (norm_nonneg _)
    exact pow_le_pow_right₀ hge (le_of_lt hi)
  have hkey : ‖α ^ n‖ = ‖∑ i ∈ Finset.range n, a i * α ^ i‖ := by
    have h : α ^ n = -(∑ i ∈ Finset.range n, a i * α ^ i) := by linear_combination hroot
    rw [h, norm_neg]
  rw [norm_pow] at hkey
  nlinarith [norm_nonneg (∑ i ∈ Finset.range n, a i * α ^ i)]

/-- **★★★★★Eisenstein の根のノルム**——`‖α‖^n = ‖a_0‖`。
`a_0` の項が一意に最大になる(他の項は `‖a_0‖·‖α‖ < ‖a_0‖`)ことによる。 -/
theorem norm_pow_eq_of_monic_root {F : Type*} [NormedField F] [IsUltrametricDist F]
    {n : ℕ} (hn : 0 < n) (a : ℕ → F) {α : F}
    (hroot : α ^ n + ∑ i ∈ Finset.range n, a i * α ^ i = 0)
    (hle : ∀ i, i < n → ‖a i‖ ≤ ‖a 0‖) (hlt : ‖a 0‖ < 1) (hne : a 0 ≠ 0) :
    ‖α‖ ^ n = ‖a 0‖ := by
  have hα1 : ‖α‖ < 1 := norm_lt_one_of_monic_root hn a hroot hle hlt hne
  have ha0pos : (0:ℝ) < ‖a 0‖ := norm_pos_iff.mpr hne
  rw [Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot hn] at hroot
  simp only [pow_zero, mul_one, zero_add] at hroot
  set S := ∑ i ∈ Finset.Ico 1 n, a i * α ^ i with hS
  have hSbound : ‖S‖ ≤ ‖a 0‖ * ‖α‖ := by
    refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
    intro i hi
    rw [Finset.mem_Ico] at hi
    rw [norm_mul, norm_pow]
    refine mul_le_mul (hle i hi.2) ?_ (by positivity) (norm_nonneg _)
    calc ‖α‖ ^ i ≤ ‖α‖ ^ 1 := pow_le_pow_of_le_one (norm_nonneg _) (le_of_lt hα1) hi.1
      _ = ‖α‖ := pow_one _
  have hSlt : ‖S‖ < ‖a 0‖ := by nlinarith [norm_nonneg α]
  have hdec : α ^ n = -(a 0) + -S := by linear_combination hroot
  have hne2 : ‖(-S : F)‖ ≠ ‖(-(a 0) : F)‖ := by
    rw [norm_neg, norm_neg]; exact ne_of_lt hSlt
  rw [← norm_pow α n, hdec, add_comm, IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm hne2,
    norm_neg, norm_neg]
  exact max_eq_right (le_of_lt hSlt)

/-! ## `ramificationIdx` を「`P^n` に入る」で下から押さえる橋

`Ideal.ramificationIdx p P = sSup {n | map p ≤ P^n}` なので、
`map p ≤ P^n` から `n ≤ e` を出すには**上に有界**であることが要る。
`Nat.sSup_of_not_bddAbove`(有界でなければ `sSup = 0`)の対偶で、
`e ≠ 0` から有界性が出る——そして `e ≠ 0` は `e·f = [K(x):K] ≥ 1` から従う。

これで「Eisenstein ⟹ 完全分岐」の道筋の (2) が通る:
`‖α‖^n = ‖π‖`(`norm_pow_eq_of_monic_root`)から `π = u·α^n`(`u` は単数)、
よって `π ∈ 𝔪_L^n`、よって `n ≤ e`、`e·f = n` と合わせて `f = 1`。 -/

/-- **`map p ≤ P^n` なら `n ≤ ramificationIdx`**(`e ≠ 0` のとき)。 -/
theorem le_ramificationIdx_of_map_le_pow {R S : Type*} [CommRing R] [CommRing S] [Algebra R S]
    {q : Ideal R} {P : Ideal S} {n : ℕ} (hne : Ideal.ramificationIdx q P ≠ 0)
    (h : Ideal.map (algebraMap R S) q ≤ P ^ n) : n ≤ Ideal.ramificationIdx q P := by
  have hbdd : BddAbove {m | Ideal.map (algebraMap R S) q ≤ P ^ m} := by
    by_contra hc
    exact hne (Nat.sSup_of_not_bddAbove hc)
  exact le_csSup hbdd h

theorem ramificationIndex_ne_zero (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    ramificationIndex K x ≠ 0 := by
  intro h
  have hmul := ramificationIndex_mul_inertiaDegree K x
  rw [h, Nat.zero_mul] at hmul
  have hpos : 0 < Module.finrank K.carrier
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := Module.finrank_pos
  omega

theorem inertiaDegree_ne_zero (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    inertiaDegree K x ≠ 0 := by
  intro h
  have hmul := ramificationIndex_mul_inertiaDegree K x
  rw [h, Nat.mul_zero] at hmul
  have hpos : 0 < Module.finrank K.carrier
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := Module.finrank_pos
  omega

/-- **★★★★★`[K(x):K] ≤ e` なら完全分岐**——`e·f = [K(x):K]` から `f ≤ 1`。
「Eisenstein ⟹ 完全分岐」の最後の一段。 -/
theorem isTotallyRamifiedAdjoin_of_finrank_le (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (h : Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≤ ramificationIndex K x) : IsTotallyRamifiedAdjoin K x := by
  have hmul := ramificationIndex_mul_inertiaDegree K x
  have hf1 : inertiaDegree K x ≤ 1 := by
    refine Nat.le_of_mul_le_mul_left ?_ (Nat.pos_of_ne_zero (ramificationIndex_ne_zero K x))
    rw [hmul, Nat.mul_one]
    exact h
  have hf0 := inertiaDegree_ne_zero K x
  show inertiaDegree K x = 1
  omega

/-! ## ★「`‖x‖^n = ‖π‖` ⟹ 完全分岐」——Eisenstein からの最後の一段 -/

/-- **★★★★★★`‖x‖^{[K(x):K]} = ‖π‖` なら `K(x)/K` は完全分岐**。

`π = u·x^n`(`u` は `‖u‖ = 1` すなわち単数)なので `π ∈ 𝔪_L^n`、よって
`n ≤ e`(`le_ramificationIdx_of_map_le_pow`)、`e·f = n` と合わせて `f = 1`。

`norm_pow_eq_of_monic_root`(Eisenstein の根のノルム)と合わせると
**「Eisenstein の根が生成する拡大は完全分岐」**が出る。 -/
theorem isTotallyRamifiedAdjoin_of_norm_pow_eq (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hnorm : ‖x‖ ^ (Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) = ‖(π : K.carrier)‖) :
    IsTotallyRamifiedAdjoin K x := by
  haveI := isLocalRing_adjoinIntegers K x
  have hnpos : 0 < Module.finrank K.carrier
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := Module.finrank_pos
  have hπle : ‖(π : K.carrier)‖ ≤ 1 := by
    have h := π.2; rw [Valued.integer.mem_iff] at h; exact h
  have hπmem : π ∈ IsLocalRing.maximalIdeal 𝒪[K.carrier] := by
    rw [hπmax]; exact Ideal.mem_span_singleton_self π
  have hπlt : ‖(π : K.carrier)‖ < 1 := by
    rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff,
      Valued.integer.isUnit_iff_norm_eq_one] at hπmem
    exact lt_of_le_of_ne hπle hπmem
  have hxlt : ‖x‖ < 1 := by
    by_contra hc
    rw [not_lt] at hc
    have h1 : (1:ℝ) ≤ ‖x‖ ^ (Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) := one_le_pow₀ hc
    rw [hnorm] at h1
    linarith
  have hxE0 : (⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ≠ 0 := by
    intro h
    have hz : ‖x‖ = 0 := by
      have hv := congrArg Subtype.val h
      show ‖x‖ = 0
      rw [show x = 0 from hv, norm_zero]
    rw [hz, zero_pow hnpos.ne'] at hnorm
    exact hπ0 (norm_eq_zero.mp hnorm.symm)
  have hxEnorm : ‖(⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ = ‖x‖ := rfl
  have hπLnorm : ‖(algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (π : K.carrier))‖
      = ‖(π : K.carrier)‖ := norm_algebraMap' _ _
  have hπnorm0 : ‖(π : K.carrier)‖ ≠ 0 := norm_ne_zero_iff.mpr hπ0
  have hu1 : ‖(algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (π : K.carrier))
      / (⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ^ (Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))‖ = 1 := by
    rw [norm_div, norm_pow, hπLnorm, hxEnorm, hnorm, div_self hπnorm0]
  have hxOle : ‖(⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ 1 := by
    rw [hxEnorm]; exact le_of_lt hxlt
  have hxOmem : (⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hxOle⟩ :
      adjoinIntegers K x) ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x) := by
    rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff, isUnit_adjoinIntegers_iff]
    show ¬ (‖(⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ = 1)
    rw [hxEnorm]; exact ne_of_lt hxlt
  have hpow0 : (⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ^ (Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ≠ 0 := pow_ne_zero _ hxE0
  have heq : (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) π
      = (⟨(algebraMap K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (π : K.carrier))
        / (⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ^ (Module.finrank K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure))), le_of_eq hu1⟩ :
        adjoinIntegers K x)
        * (⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hxOle⟩ :
          adjoinIntegers K x) ^ (Module.finrank K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) := by
    apply Subtype.ext
    show algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (π : K.carrier) = _
    push_cast
    rw [div_mul_cancel₀ _ hpow0]
  have hmem : (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) π
      ∈ (IsLocalRing.maximalIdeal (adjoinIntegers K x)) ^ (Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) := by
    rw [heq]
    exact Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hxOmem _)
  have hmap : Ideal.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      (IsLocalRing.maximalIdeal 𝒪[K.carrier])
      ≤ (IsLocalRing.maximalIdeal (adjoinIntegers K x)) ^ (Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) := by
    rw [hπmax, Ideal.map_span, Ideal.span_le, Set.image_singleton, Set.singleton_subset_iff]
    exact hmem
  exact isTotallyRamifiedAdjoin_of_finrank_le K x
    (le_ramificationIdx_of_map_le_pow (ramificationIndex_ne_zero K x) hmap)

end ABC3.Found.PGC
