import ABC3.Found.PGC.AbsClosureModules
import ABC3.Found.PGC.GaloisTransferContinuous
import ABC3.Found.PGC.ResidueCardinality

/-!
# [pGC] Proposition 2.2 —— 自由なデータ引数版は偽。実物に固定した形とその還元

原文 (pGC 物理 p.5, Proposition 2.2):

> Suppose that we are given the following group-theoretic data: the topological group Γ_K,
> together with the indexed filtration Γ_K^v for all v > 0. Then the Γ_K-modules O[scr]_K[bar], and
> K[bar]∧ can be recovered group-theoretically from this group-theoretic data.

同, Proof(`0_Source` の `.txt` を切り出したもの。抽出器の字形の潰れは直していない):

> To obtain OK, we simply apply the above discussion to finite extensions L of K,
> and use Proposition 2.1. Note that when we pass to an open subgroup ΓL ⊆ΓK, since the
> upper numbering is not compatible with passage to subgroups, one must first convert to
> the lower numbering (which is compatible with passage to subgroups), and then convert
> back to the upper numbering for ΓL.

## ★★★1. `Skeleton/PGC/Section2.lean::prop_2_2` の形は**偽である**

スケルトンの形は `𝒪_{K̄}`・`K̄^∧` を**自由な型族**として取る:

```lean
theorem prop_2_2 (_RF : RamificationFiltration p)
    (IntKbar CompKbar : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (IntKbar K)] [∀ K, DistribMulAction K.absGal (IntKbar K)]
    [∀ K, AddCommGroup (CompKbar K)] [∀ K, DistribMulAction K.absGal (CompKbar K)] :
    RecoverableAsAddModule IntKbar ∧ RecoverableAsAddModule CompKbar
```

★これは **D13**(`Skeleton/PGC/Section1.lean` の Prop 1.2 の docstring)が
`∀ RD : ResidueCardinality p` について記録したのと**同じ退化**である。
`Check/PGC/Prop22Degenerate.lean`(2026-09-05)は作用を `SMul` から
`DistribMulAction` に強めて 1 つ目の反例を塞いだが、★**それでもまだ偽**である。

理由は病的な作用ではなく、**型族が `PAdicLocalField p` の「項」の関数**であること:

* `Γ_K ≅ Γ_{K'}`(位相群として)なのに `K ≠ K'` である項の対が実在する
  (`Check/PGC/Prop12Degenerate.lean` の `twistedField p` と `selfField p`。
  台の型を `ℚ_[p]` のままにして体構造だけを `x ↦ -x` で捻ったもの)。
* そこで「`K = twistedField p` のときだけ `ℤ`、他は `0`」という型族を取れば
  `Obj K ≃+ Obj K'` は `ℤ ≃+ 0` を要求して落ちる。作用は**自明な作用**でよいので
  `DistribMulAction` の公理は全部満たされている。

本ファイル §1 がその**抽象核**(`condSubgroup`・`false_of_addEquiv_condSubgroup`)を置き、
具体層(反例そのもの)は `Check/PGC/Prop22FreeForm.lean::prop_2_2_free_form_false` に在る。

★★**落とした条件は「同型不変性」である**——D13 と同じ。原典は `𝒪_{K̄}`・`K̄^∧` という
**実物**を語っているのであって、任意の型族を語ってはいない。

## 2. 実物に固定した形

実物は `Found/PGC/AbsClosureModules.lean` に構成済みである:

```
IntKbar K  = ↥(absClosureInt K)   = 𝒪_{K̄}
CompKbar K = closureCompletion K  = K̄^∧ = ℂ_K
```

本ファイルはこれらについて `IntKbarRecoverable` / `CompKbarRecoverable` を
**ただ 1 つの仮説**に還元する(還元自体は `sorry` 無し):

```
IsometricallyRecoverableClosure p :=
  ∀ {K K'} (α : Γ_K ≃ₜ* Γ_{K'}), ∃ φ : K̄ ≃+ K̄', (φ は等長) ∧ (φ は α-同変)
```

すなわち ★**「Proposition 2.1 の同型を等長に取れる」**。
分岐フィルトレーション `Γ_K^v` が担っているのは**まさにこの一点**である
——原文が上付き/下付き番号付けの変換を経由して `𝒪_L ⊆ L` を群論的に切り出すのは、
`K̄` の付値(= ノルム)を群論的データから復元するためだからである。

★★**Proposition 2.1(`Found/PGC/CountableGenerators.lean::prop_2_1`、`sorry` 無し)は
「同型が在る」までしか言わない。**`recoverableAsAddModule_closure_of_isometric` は
本仮説がその **真の強化**であることを示す(仮説 ⇒ Prop 2.1)。

## 3. ★原典より短い道

原文の証明は「有限次拡大 L/K へ降りて Prop 2.1 を使い、上付き→下付き→上付きと
番号付けを往復する」というものである。本ファイルはそこを**一切通らない**:

* 有限次拡大への降下も、Herbrand の変換も使わない。
* 代わりに ★**「等長な α-同変加法同型は単位球へ制限され、かつ完備化へ延びる」**
  という 2 本の抽象核(§2)だけで `𝒪_{K̄}` と `ℂ_K` の**両方**を同時に得る。
  原文が `K̄^∧` について別途「`𝒪_{K̄}` の p 進完備化を `ℤ_p` 上 `ℚ_p` と
  テンソルする」と述べる段は、`addEquivCompletion` 1 本に吸収される。

★これは逸脱ではない:得られる主張は同一で、証明の道筋だけが違う。

## 4. ★非空虚性(仮説を空虚に満たしていないことの検査)

`IsometricTransport α`(α ごとの形)は、次の 2 つの場合に**無条件に成り立つ**:

| 場合 | 補題 | 使う φ |
|---|---|---|
| α が内部自己同型 `g ↦ c g c⁻¹` | `isometricTransport_inner` | `x ↦ c x` |
| α が体の同型 `β : K ≃ₐ[ℚ_p] K'` から来る | `isometricTransport_galContinuousMulEquiv` | `extendToClosure β` |

★★**2 番目は、上の反例が使うのとまったく同じ α である。**
すなわち ★**自由版が落ちるその α において、実物版は成り立つ**
(`Check/PGC/Prop22FreeForm.lean::intKbar_transport_twisted`)。
仮説は空虚でも、反例の場所を避けてもいない。

そのために新しく証明したのが

* `norm_extendToClosure` —— ★**ℚ_p-代数同型の代数閉包への延長はスペクトルノルムを保つ**。

である。既存の `norm_algEquiv_carrier`(`Found/PGC/ResidueCardinality.lean`、
有限次の底の間)を `spectralNorm_unique_field_norm_ext` で代数閉包へ持ち上げた。
★木にも mathlib にも無かった(実測:
`grep -n "ABC3.Found.PGC" .cache/decl-index.txt | grep -i norm | grep -i "algEquiv\|equiv"`
は `norm_algEquiv_carrier`(底のみ)しか返さない)。

## 5. 逸脱の記録

* `Skeleton/PGC/Section2.lean` は**書き換えていない**(D13 と同じく、スケルトンの
  修正は本体が別の波でやる)。本ファイルは Found 側に「実物に固定した形」を出すだけである。
* `prop_2_2` 本体を**無条件には証明していない**。残るのは
  `IsometricallyRecoverableClosure`(= 分岐フィルトレーションが担う一点)である。
  ★これは「数学が足りない」側の穴であって、配管の穴ではない。
* 原文の `K̄^∧` を「`𝒪_{K̄}` の p 進完備化 ⊗_{ℤ_p} ℚ_p」ではなく
  `K^al` のノルム完備化として扱う点は `AbsClosureModules.lean` の同定
  (`mem_pow_p_mul_absClosureInt`)をそのまま踏襲する。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open UniformSpace

/-! ## §1 抽象核 A —— 「自由な型族」版が偽であることの核

★ここには分岐も付値も Galois も `PAdicLocalField` も出てこない。
出てくるのは「命題 1 つ」と「`ℤ` の部分群」だけである。 -/

section FreeFamilyCore

/-- 命題 `P` が真なら `ℤ`(`⊤`)、偽なら `0`(`⊥`)になる `ℤ` の部分群。

★`if` を使わないので `DecidableEq` も `Classical` も要らない。 -/
def condSubgroup (P : Prop) : AddSubgroup ℤ where
  carrier := {n : ℤ | ¬ P → n = 0}
  zero_mem' _ := rfl
  add_mem' ha hb hP := by rw [ha hP, hb hP, add_zero]
  neg_mem' ha hP := by rw [ha hP, neg_zero]

theorem nontrivial_condSubgroup {P : Prop} (hP : P) : Nontrivial ↥(condSubgroup P) :=
  ⟨⟨⟨0, fun h => absurd hP h⟩, ⟨1, fun h => absurd hP h⟩,
    fun h => by simpa using congrArg Subtype.val h⟩⟩

theorem subsingleton_condSubgroup {P : Prop} (hP : ¬ P) : Subsingleton ↥(condSubgroup P) :=
  ⟨fun x y => Subtype.ext ((x.2 hP).trans (y.2 hP).symm)⟩

/-- 自明な作用。★`DistribMulAction` の公理は**すべて満たす**——
`Check/PGC/Prop22Degenerate.lean` が塞いだ「公理ゼロの `SMul`」型の病理では**ない**。 -/
@[implicit_reducible] def constDistribMulAction (M A : Type*) [Monoid M] [AddCommGroup A] :
    DistribMulAction M A where
  smul _ a := a
  one_smul _ := rfl
  mul_smul _ _ _ := rfl
  smul_zero _ := rfl
  smul_add _ _ _ := rfl

/-- **★★★抽象核 A** —— 「真の点」と「偽の点」の間に加法同型は無い。

これが「自由な型族 `Obj : ι → Type` を全称量化した回復可能性」を壊す核である:
`ι` の 2 点 `i ≠ j` を区別する族 `Obj k := ↥(condSubgroup (k = i))` を取れば
`Obj i ≃+ Obj j` は `ℤ ≃+ 0` を要求する。 -/
theorem false_of_addEquiv_condSubgroup {P Q : Prop} (hP : P) (hQ : ¬ Q)
    (φ : ↥(condSubgroup P) ≃+ ↥(condSubgroup Q)) : False := by
  haveI := subsingleton_condSubgroup hQ
  obtain ⟨a, b, hab⟩ := (nontrivial_condSubgroup hP).exists_pair_ne
  exact hab (φ.injective (Subsingleton.elim _ _))

end FreeFamilyCore

/-! ## §2 抽象核 B —— 等長な同変加法同型は単位球へ制限され、完備化へ延びる

★ここにも分岐・付値・Galois・`PAdicLocalField` は 1 語も出てこない。
出てくるのは「ノルム付き加法群」と「モノイド作用」だけである。 -/

section RestrictExtendCore

variable {A B : Type*}

/-- **抽象核 B1** —— 加法同型が 2 つの部分群を対応させるなら、部分群の間の加法同型に制限される。 -/
def addEquivOfMemIff [AddCommGroup A] [AddCommGroup B]
    {SA SB : Type*} [SetLike SA A] [AddSubgroupClass SA A] [SetLike SB B] [AddSubgroupClass SB B]
    (φ : A ≃+ B) (sa : SA) (sb : SB) (h : ∀ x : A, x ∈ sa ↔ φ x ∈ sb) : ↥sa ≃+ ↥sb where
  toFun x := ⟨φ x, (h x).mp x.2⟩
  invFun y := ⟨φ.symm y, (h _).mpr (by rw [AddEquiv.apply_symm_apply]; exact y.2)⟩
  left_inv x := Subtype.ext (φ.symm_apply_apply _)
  right_inv y := Subtype.ext (φ.apply_symm_apply _)
  map_add' x y := Subtype.ext (φ.map_add _ _)

@[simp] theorem coe_addEquivOfMemIff [AddCommGroup A] [AddCommGroup B]
    {SA SB : Type*} [SetLike SA A] [AddSubgroupClass SA A] [SetLike SB B] [AddSubgroupClass SB B]
    (φ : A ≃+ B) (sa : SA) (sb : SB) (h : ∀ x : A, x ∈ sa ↔ φ x ∈ sb) (x : ↥sa) :
    ((addEquivOfMemIff φ sa sb h x : ↥sb) : B) = φ (x : A) := rfl

/-- **抽象核 B1'** —— ノルムを保つ加法同型は「ノルム ≤ 1 の元」の集合を対応させる。 -/
theorem memIff_of_norm_eq [SeminormedAddCommGroup A] [SeminormedAddCommGroup B]
    {SA SB : Type*} [SetLike SA A] [SetLike SB B]
    (φ : A ≃+ B) (hφ : ∀ x, ‖φ x‖ = ‖x‖) (sa : SA) (sb : SB)
    (hsa : ∀ x : A, x ∈ sa ↔ ‖x‖ ≤ 1) (hsb : ∀ y : B, y ∈ sb ↔ ‖y‖ ≤ 1) (x : A) :
    x ∈ sa ↔ φ x ∈ sb := by rw [hsa, hsb, hφ]

/-- **抽象核 B1''** —— 制限した同型も同変。 -/
theorem addEquivOfMemIff_smul {M M' : Type*} [Monoid M] [Monoid M'] [AddCommGroup A]
    [AddCommGroup B] [DistribMulAction M A] [DistribMulAction M' B]
    {SA SB : Type*} [SetLike SA A] [AddSubgroupClass SA A] [SetLike SB B] [AddSubgroupClass SB B]
    (φ : A ≃+ B) (sa : SA) (sb : SB) (h : ∀ x : A, x ∈ sa ↔ φ x ∈ sb)
    [SMul M ↥sa] [SMul M' ↥sb] (α : M → M')
    (hA : ∀ (g : M) (x : ↥sa), ((g • x : ↥sa) : A) = g • (x : A))
    (hB : ∀ (g' : M') (y : ↥sb), ((g' • y : ↥sb) : B) = g' • (y : B))
    (heq : ∀ (g : M) (x : A), φ (g • x) = α g • φ x) (g : M) (x : ↥sa) :
    addEquivOfMemIff φ sa sb h (g • x) = α g • addEquivOfMemIff φ sa sb h x :=
  Subtype.ext (by rw [coe_addEquivOfMemIff, hA, hB, coe_addEquivOfMemIff, heq])

/-- **抽象核 B2** —— 両向き連続な加法同型は完備化の加法同型へ延びる。

mathlib には環同型の版(`UniformSpace.Completion.mapRingEquiv`)しか無い。
加法同型の版は `AddMonoidHom.completion` を両向きに使って組む。 -/
noncomputable def addEquivCompletion [AddCommGroup A] [UniformSpace A] [IsUniformAddGroup A]
    [AddCommGroup B] [UniformSpace B] [IsUniformAddGroup B]
    (φ : A ≃+ B) (hc : Continuous φ) (hc' : Continuous φ.symm) :
    Completion A ≃+ Completion B where
  toFun := AddMonoidHom.completion (φ : A →+ B) hc
  invFun := AddMonoidHom.completion (φ.symm : B →+ A) hc'
  left_inv z := by
    refine Completion.induction_on z (isClosed_eq ?_ continuous_id) ?_
    · exact (AddMonoidHom.continuous_completion _ hc').comp
        (AddMonoidHom.continuous_completion _ hc)
    · intro a
      rw [AddMonoidHom.completion_coe, AddMonoidHom.completion_coe]; simp
  right_inv z := by
    refine Completion.induction_on z (isClosed_eq ?_ continuous_id) ?_
    · exact (AddMonoidHom.continuous_completion _ hc).comp
        (AddMonoidHom.continuous_completion _ hc')
    · intro b
      rw [AddMonoidHom.completion_coe, AddMonoidHom.completion_coe]; simp
  map_add' := map_add _

@[simp] theorem addEquivCompletion_coe [AddCommGroup A] [UniformSpace A] [IsUniformAddGroup A]
    [AddCommGroup B] [UniformSpace B] [IsUniformAddGroup B]
    (φ : A ≃+ B) (hc : Continuous φ) (hc' : Continuous φ.symm) (a : A) :
    addEquivCompletion φ hc hc' (a : Completion A) = ((φ a : B) : Completion B) :=
  AddMonoidHom.completion_coe _ hc a

/-- **抽象核 B2'** —— 完備化へ延ばした同型も同変。稠密性と両辺の連続性だけで出る。 -/
theorem addEquivCompletion_smul {M M' : Type*} [Monoid M] [Monoid M']
    [AddCommGroup A] [UniformSpace A] [IsUniformAddGroup A]
    [AddCommGroup B] [UniformSpace B] [IsUniformAddGroup B]
    [DistribMulAction M A] [UniformContinuousConstSMul M A]
    [DistribMulAction M' B] [UniformContinuousConstSMul M' B]
    (φ : A ≃+ B) (hc : Continuous φ) (hc' : Continuous φ.symm) (α : M → M')
    (heq : ∀ (g : M) (x : A), φ (g • x) = α g • φ x) (g : M) (z : Completion A) :
    addEquivCompletion φ hc hc' (g • z) = α g • addEquivCompletion φ hc hc' z := by
  refine Completion.induction_on z (isClosed_eq ?_ ?_) ?_
  · exact (AddMonoidHom.continuous_completion _ hc).comp (continuous_const_smul g)
  · exact (continuous_const_smul (α g)).comp (AddMonoidHom.continuous_completion _ hc)
  · intro a
    rw [← Completion.coe_smul, addEquivCompletion_coe, addEquivCompletion_coe, heq,
      Completion.coe_smul]

end RestrictExtendCore

/-! ## §3 具体層 —— 実物 `𝒪_{K̄}`・`ℂ_K` への代入 -/

variable {p : ℕ} [Fact p.Prime]

/-- **α ごとの「等長版 Proposition 2.1」**。

`α : Γ_K ≃ₜ* Γ_{K'}` に沿う `K̄ ≃+ K̄'` の同変同型で、しかも**ノルムを保つ**ものの存在。
★分岐フィルトレーション `Γ_K^v` が担う内容はここに集約されている。 -/
def IsometricTransport {K K' : PAdicLocalField p}
    (α : ContinuousMulEquiv K.absGal K'.absGal) : Prop :=
  ∃ φ : K.closure ≃+ K'.closure, (∀ x, ‖φ x‖ = ‖x‖) ∧
    ∀ (g : K.absGal) (x : K.closure), φ (g • x) = (α.toMulEquiv g) • φ x

/-- **等長版 Proposition 2.1**(すべての `α` について)。 -/
def IsometricallyRecoverableClosure (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : ContinuousMulEquiv K.absGal K'.absGal), IsometricTransport α

/-- ★本仮説は **Proposition 2.1 の真の強化**である(等長性を落とせば Prop 2.1)。 -/
theorem recoverableAsAddModule_closure_of_isometric (h : IsometricallyRecoverableClosure p) :
    RecoverableAsAddModule (p := p) (fun K => K.closure) := by
  intro K K' α
  obtain ⟨φ, -, hequiv⟩ := h α
  exact ⟨φ, hequiv⟩

/-! ### `𝒪_{K̄}` —— 単位球への制限(抽象核 B1 の代入) -/

/-- **α ごとの `𝒪_{K̄}` の移送**。 -/
theorem intKbar_transport_of_isometricTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IsometricTransport α) :
    ∃ φ : IntKbar K ≃+ IntKbar K',
      ∀ (g : K.absGal) (x : IntKbar K), φ (g • x) = (α.toMulEquiv g) • φ x := by
  obtain ⟨ψ, hnorm, hequiv⟩ := h
  have hmem : ∀ x : K.closure, x ∈ absClosureInt K ↔ ψ x ∈ absClosureInt K' := by
    intro x; rw [mem_absClosureInt, mem_absClosureInt, hnorm]
  refine ⟨addEquivOfMemIff ψ (absClosureInt K) (absClosureInt K') hmem, fun g x => ?_⟩
  exact Subtype.ext (by
    rw [coe_addEquivOfMemIff, coe_smul_absClosureInt, coe_smul_absClosureInt,
      coe_addEquivOfMemIff, ← smul_closure_def, hequiv, smul_closure_def])

/-- **`𝒪_{K̄}` は群論的に回復できる**——等長版 Prop 2.1 を仮定すれば。 -/
theorem intKbar_recoverable_of_isometric (h : IsometricallyRecoverableClosure p) :
    IntKbarRecoverable (p := p) := fun α => intKbar_transport_of_isometricTransport (h α)

/-! ### `ℂ_K = K̄^∧` —— 完備化への延長(抽象核 B2 の代入) -/

/-- **α ごとの `ℂ_K` の移送**。 -/
theorem compKbar_transport_of_isometricTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IsometricTransport α) :
    ∃ φ : CompKbar K ≃+ CompKbar K',
      ∀ (g : K.absGal) (x : CompKbar K), φ (g • x) = (α.toMulEquiv g) • φ x := by
  obtain ⟨ψ, hnorm, hequiv⟩ := h
  have hsymm : ∀ y : K'.closure, ‖ψ.symm y‖ = ‖y‖ := by
    intro y; rw [← hnorm (ψ.symm y), ψ.apply_symm_apply]
  have hc : Continuous ψ := (AddMonoidHomClass.isometry_of_norm ψ hnorm).continuous
  have hc' : Continuous ψ.symm := (AddMonoidHomClass.isometry_of_norm ψ.symm hsymm).continuous
  exact ⟨addEquivCompletion ψ hc hc',
    addEquivCompletion_smul ψ hc hc' (fun g => α.toMulEquiv g) hequiv⟩

/-- **`ℂ_K` は群論的に回復できる**——等長版 Prop 2.1 を仮定すれば。 -/
theorem compKbar_recoverable_of_isometric (h : IsometricallyRecoverableClosure p) :
    CompKbarRecoverable (p := p) := fun α => compKbar_transport_of_isometricTransport (h α)

/-- ★★★**[pGC] Proposition 2.2(実物に固定した形)** ——
`IsometricallyRecoverableClosure`(= 分岐フィルトレーションが担う一点)への還元。

★自由な型族版(`Skeleton/PGC/Section2.lean::prop_2_2` の形)は**偽**である
(`Check/PGC/Prop22FreeForm.lean::prop_2_2_free_form_false`)。 -/
theorem prop_2_2_real_of_isometric (h : IsometricallyRecoverableClosure p) :
    IntKbarRecoverable (p := p) ∧ CompKbarRecoverable (p := p) :=
  ⟨intKbar_recoverable_of_isometric h, compKbar_recoverable_of_isometric h⟩

def prop_2_2_real_of_isometric.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-! ## §4 非空虚性 —— 仮説が実際に満たされる場合 -/

/-- `Γ_K` の内部自己同型(位相群の同型)。 -/
noncomputable def innerAbsGalEquiv (K : PAdicLocalField p) (c : K.absGal) :
    ContinuousMulEquiv K.absGal K.absGal where
  toFun g := c * g * c⁻¹
  invFun g := c⁻¹ * g * c
  left_inv g := by group
  right_inv g := by group
  map_mul' a b := by group
  continuous_toFun := by fun_prop
  continuous_invFun := by fun_prop

@[simp] theorem innerAbsGalEquiv_apply (K : PAdicLocalField p) (c g : K.absGal) :
    innerAbsGalEquiv K c g = c * g * c⁻¹ := rfl

theorem isometricTransport_refl (K : PAdicLocalField p) :
    IsometricTransport (ContinuousMulEquiv.refl K.absGal) :=
  ⟨AddEquiv.refl _, fun _ => rfl, fun _ _ => rfl⟩

/-- ★**内部自己同型に対しては仮説が無条件に成り立つ**(φ は `c` の作用そのもの)。 -/
theorem isometricTransport_inner (K : PAdicLocalField p) (c : K.absGal) :
    IsometricTransport (innerAbsGalEquiv K c) := by
  refine ⟨(c : K.closure ≃ₐ[K.carrier] K.closure).toRingEquiv.toAddEquiv, ?_, ?_⟩
  · intro x; exact norm_smul_closure K c x
  · intro g x
    show c (g x) = (c * g * c⁻¹) • (c x)
    show c (g x) = (c * g * c⁻¹) (c x)
    show c (g x) = c (g (c⁻¹ (c x)))
    congr 2
    exact (c.symm_apply_apply x).symm

/-- ★★**`ℚ_p`-代数同型の代数閉包への延長はスペクトルノルムを保つ**(逆向き)。

底の間の等長性 `norm_algEquiv_carrier`(`Found/PGC/ResidueCardinality.lean`)を
`spectralNorm_unique_field_norm_ext` で代数閉包へ持ち上げたもの。
★木にも mathlib にも無かった。 -/
theorem norm_extendToClosure_symm {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (y : K'.closure) :
    ‖(extendToClosure β).symm y‖ = ‖y‖ := by
  have hinj : Function.Injective
      (((extendToClosure β).symm : K'.closure ≃+* K.closure) : K'.closure →+* K.closure) :=
    (extendToClosure β).symm.injective
  set f : AbsoluteValue K'.closure ℝ :=
    (NormedField.toAbsoluteValue K.closure).comp hinj with hf
  have hext : ∀ a : K'.carrier, f (algebraMap K'.carrier K'.closure a) = ‖a‖ := by
    intro a
    show ‖(extendToClosure β).symm (algebraMap K'.carrier K'.closure a)‖ = ‖a‖
    have hb : (extendToClosure β) (algebraMap K.carrier K.closure (β.symm a))
        = algebraMap K'.carrier K'.closure a := by
      rw [extendToClosure_algebraMap, β.apply_symm_apply]
    rw [← hb, RingEquiv.symm_apply_apply, norm_algebraMap_closure, norm_algEquiv_carrier β.symm]
  have h := spectralNorm_unique_field_norm_ext (K := K'.carrier) (L := K'.closure) (f := f) hext y
  rw [← norm_eq_spectralNorm_closure] at h
  exact h

/-- ★★同上(順向き)。 -/
theorem norm_extendToClosure {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (x : K.closure) :
    ‖extendToClosure β x‖ = ‖x‖ := by
  have h := norm_extendToClosure_symm β (extendToClosure β x)
  rw [RingEquiv.symm_apply_apply] at h
  exact h.symm

/-- ★★★**体の同型から来る `α` については仮説が無条件に成り立つ**。

★これは `Check/PGC/Prop12Degenerate.lean` の `twistedGalEquiv`——
すなわち**自由な型族版の反例が使う α そのもの**を含む。
自由版が落ちるその場所で、実物版は成り立つ。 -/
theorem isometricTransport_galContinuousMulEquiv {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    IsometricTransport (galContinuousMulEquiv β) := by
  refine ⟨(extendToClosure β : K.closure ≃+* K'.closure).toAddEquiv, norm_extendToClosure β, ?_⟩
  intro g x
  show extendToClosure β (g x)
    = (conjGalOfEquiv β (extendToClosure β) (extendToClosure_algebraMap β) g)
        (extendToClosure β x)
  rw [conjGalOfEquiv_apply, RingEquiv.symm_apply_apply]
  rfl

/-- ★体の同型から来る `α` については `𝒪_{K̄}` の移送が**無条件に**存在する。 -/
theorem intKbar_transport_galContinuousMulEquiv {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    ∃ φ : IntKbar K ≃+ IntKbar K', ∀ (g : K.absGal) (x : IntKbar K),
      φ (g • x) = ((galContinuousMulEquiv β).toMulEquiv g) • φ x :=
  intKbar_transport_of_isometricTransport (isometricTransport_galContinuousMulEquiv β)

/-- ★体の同型から来る `α` については `ℂ_K` の移送が**無条件に**存在する。 -/
theorem compKbar_transport_galContinuousMulEquiv {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    ∃ φ : CompKbar K ≃+ CompKbar K', ∀ (g : K.absGal) (x : CompKbar K),
      φ (g • x) = ((galContinuousMulEquiv β).toMulEquiv g) • φ x :=
  compKbar_transport_of_isometricTransport (isometricTransport_galContinuousMulEquiv β)

def intKbar_transport_galContinuousMulEquiv.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

end ABC3.Found.PGC
