import ABC3.Found.FrdI.Prop16
import Mathlib.CategoryTheory.Action

/-!
# `Aut-ample` は `Definition 1.3` から出るか —— `[FrdI] Proposition 1.6, (v)` の `⟸`

★`Proposition 1.6, (v)` の `base-trivial` / `metrically trivial` の `⟸` には
「**底を指定した同型を取り直せる**」ことが要る(`Gap/FrdI/Section1.lean` の `Gap_1_6_v`)。
それは `Definition 1.2, (iv)` の `Aut-ample` そのものである。

## ★★反例の設計

`𝒟` を **1 対象の亜群** `SingleObj Γ` とし、`𝒞` を **`Γ`-集合 `S` 上に広げた `𝔽_Φ`**
——すなわち作用亜群 `Γ ⋉ S` 上の `𝔽_{Φ∘π}` を、`π : Γ ⋉ S ⥤ SingleObj Γ` を通して
**`𝒟` 上の pre-Frobenioid と見たもの**——とする。すると:

- ★`𝒟` は 1 対象なので **`BaseIsomorphic` は常に真**。したがって
  `base-trivial` ⟺ **`𝒞` の対象がすべて同型** ⟺ ★**`Γ` の作用が推移的**
- ★`Aut_𝒞(s) → Aut_𝒟(*) = Γ` の像は **`Stab(s)`**。したがって
  ★★**`Aut-ample` は `Stab(s) = Γ` と同値**
- ★`Definition 1.3, (i), (b)`(`preStepSpan`)は
  ★★**`Γ = ⋃_z z·Stab(s)·z⁻¹`** を要求する

★★**したがって「推移的で、各元が不動点を持ち、固定部分群が真部分群」な作用があればよい。**
★**有限台の置換のなす群 `Γ` の `ℕ` への作用**がそれである ——
台が有限なので必ず不動点があり、互換で推移的で、`Stab(0) ≠ Γ`。

## ★測定 —— ここまでで確定したこと(2026-08-16)

★**確定した(このファイルで証明済み)**:
- `FinSuppPerm` が部分群であること
- ★**作用が推移的**(`gam_transitive`)
- ★**どの元も不動点を持つ**(`gam_has_fixedPoint`)
- ★**固定部分群が真部分群**(`gam_stab_proper`)
- ★★★**`Γ = ⋃_z z·Stab(0)·z⁻¹`**(`gam_conj_into_stab`)——
  ★**有限群では起こり得ない**(Jordan)ので、無限群が要る理由もここにある
- 2 つの底圏が connected・totally epimorphic

★**未確定(残る仕事)**:
- ★★**`𝒞 = 𝔽_{Φ∘π}` を `𝒟 = SingleObj Γ` 上の pre-Frobenioid と見たとき、
  `Definition 1.3` の全条件を満たすか。**
  ★`Base` を含まない条項は `Dact` 上の版(`elemFrob_frobenioidCore`)と
  **同じ主張になる**(`𝒟`・`Dact` はどちらも亜群なので `IsBaseIsomorphism` は常に真)。
  ★**作り直しが要るのは `baseSurj`・`preStepSpan`・`plBkEquiv` の 3 条**と見ている
  (`IsPullBack` は紙の上では両者で一致する —— `𝒟` が 1 対象なので
  ファイバー積の `h` 成分が `g` から一意に決まり、条件が「`deg = 1` かつ `Div = 0`」に落ちる)。
- ★その上で `𝒞' = 𝒞 ×_𝒟 (終圏)` を取り、
  ★**`A` は base-trivial だが `A'` はそうでない**ことを示す。

★★**`sorry` は置かない。** 未確定のものは**書かない**。
-/

namespace ABC3.Check.FrdI

open CategoryTheory ABC3.Found.FrdI

universe v u w v3 u3

/-! ## ★有限台の置換のなす群 -/

/-- ★**有限台の置換**のなす `Equiv.Perm ℕ` の部分群。 -/
def FinSuppPerm : Subgroup (Equiv.Perm ℕ) where
  carrier := {σ | {n | σ n ≠ n}.Finite}
  mul_mem' {a b} ha hb := by
    refine Set.Finite.subset (hb.union ha) ?_
    intro n hn
    by_contra hcon
    simp only [Set.mem_union, Set.mem_setOf_eq, not_or, not_not] at hcon
    exact hn (by simp only [Equiv.Perm.mul_apply, hcon.1, hcon.2])
  one_mem' := by simp
  inv_mem' {a} ha := by
    refine Set.Finite.subset ha ?_
    intro n hn
    simp only [Set.mem_setOf_eq] at hn ⊢
    intro h
    refine hn ?_
    have hc := congrArg (fun x : ℕ => (a⁻¹ : Equiv.Perm ℕ) x) h
    simpa using hc.symm

/-- ★反例に使う群 —— 有限台の置換のなす群。 -/
abbrev Gam : Type := FinSuppPerm

instance : MulAction Gam ℕ := inferInstance

@[simp] theorem gam_smul (γ : Gam) (n : ℕ) : γ • n = (γ : Equiv.Perm ℕ) n := rfl

/-- ★★**作用は推移的** —— 互換は有限台。 -/
theorem gam_transitive (a b : ℕ) : ∃ γ : Gam, γ • a = b := by
  refine ⟨⟨Equiv.swap a b, ?_⟩, ?_⟩
  · refine Set.Finite.subset (Set.toFinite ({a, b} : Set ℕ)) ?_
    intro n hn
    simp only [Set.mem_setOf_eq] at hn
    by_contra hcon
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff, not_or] at hcon
    exact hn (Equiv.swap_apply_of_ne_of_ne hcon.1 hcon.2)
  · show Equiv.swap a b a = b
    simp

/-- ★★**どの元も不動点を持つ** —— 台が有限で `ℕ` が無限だから。 -/
theorem gam_has_fixedPoint (γ : Gam) : ∃ n : ℕ, γ • n = n := by
  by_contra hcon
  have hsub : (Set.univ : Set ℕ) ⊆ {n | (γ : Equiv.Perm ℕ) n ≠ n} :=
    fun n _ h => hcon ⟨n, h⟩
  exact Set.infinite_univ (Set.Finite.subset γ.2 hsub)

/-- ★**固定部分群は真部分群** —— `0` を動かす有限台の置換がある。 -/
theorem gam_stab_proper : ∃ γ : Gam, γ • 0 ≠ 0 := by
  obtain ⟨γ, hγ⟩ := gam_transitive 0 1
  exact ⟨γ, by rw [hγ]; exact Nat.one_ne_zero⟩

/-- ★★★**`Γ` は固定部分群の共役の合併である** ——
`Γ = ⋃_z z·Stab(0)·z⁻¹`。

★★**これが `Definition 1.3, (i), (b)`(`preStepSpan`)を通す鍵**である。
有限群ではこれは真部分群では起こり得ない(Jordan)ので、
★**無限群でなければならない。** -/
theorem gam_conj_into_stab (α : Gam) : ∃ z : Gam, (z⁻¹ * α * z) • (0 : ℕ) = 0 := by
  obtain ⟨n, hn⟩ := gam_has_fixedPoint α
  obtain ⟨z, hz⟩ := gam_transitive 0 n
  refine ⟨z, ?_⟩
  rw [mul_smul, mul_smul, hz, hn, ← hz, inv_smul_smul]

/-! ## ★2 つの底圏

`𝒟 = SingleObj Γ`(1 対象)と、作用亜群 `Γ ⋉ ℕ`。
★**`𝒞` は後者の上の `𝔽`、しかし pre-Frobenioid としては前者の上に置く。** -/

/-- ★1 対象の底圏。 -/
abbrev DSingle : Type := SingleObj Gam

/-- ★作用亜群 `Γ ⋉ ℕ`。 -/
abbrev DAct : Type := ActionCategory Gam ℕ

/-- ★亜群は totally epimorphic(すべての射が同型)。 -/
theorem isTotallyEpimorphic_of_groupoid (X : Type*) [Groupoid X] :
    IsTotallyEpimorphic X := fun _ _ _ => inferInstance

theorem dsingle_totEpi : IsTotallyEpimorphic DSingle :=
  isTotallyEpimorphic_of_groupoid _

theorem dact_totEpi : IsTotallyEpimorphic DAct :=
  isTotallyEpimorphic_of_groupoid _

instance dsingle_isConnected : IsConnected DSingle := by
  refine IsConnected.of_induct (j₀ := (SingleObj.star Gam)) ?_
  intro p hp0 _ A
  obtain rfl : A = SingleObj.star Gam := Subsingleton.elim _ _
  exact hp0

/-- ★★**作用亜群は connected** —— 作用が推移的だから。 -/
instance dact_isConnected : IsConnected DAct := by
  refine IsConnected.of_induct (j₀ := (⟨(), (0 : ℕ)⟩ : DAct)) ?_
  intro p hp0 hstep A
  obtain ⟨γ, hγ⟩ := gam_transitive 0 A.back
  have hobj : A = (⟨(), γ • (0 : ℕ)⟩ : DAct) := by
    obtain ⟨⟨⟩, a⟩ := A
    simp only [ActionCategory.back] at hγ
    rw [hγ]
  rw [hobj]
  exact (hstep (ActionCategory.homOfPair (γ • (0 : ℕ)) γ)).mp
    (by simpa using hp0)

/-! ## ★`𝒟` 上の pre-Frobenioid としての `𝒞 = 𝔽_{Φ∘π}` -/

/-- ★底圏の射はすべて同型なので FSM 射。 -/
theorem dact_fsm {A B : DAct} (α : B ⟶ A) (_ : IsFSMMorphism α) :
    IsFSMMorphism ((ActionCategory.π Gam ℕ).map α) :=
  isFSMMorphism_of_isIso _

/-- ★1 対象の底圏上の `Φ = ℕ`(定数)。 -/
abbrev AgPhi : MonoidOn.{0, 0, 0} DSingle := MonoidOn.const DSingle ℕ

/-- ★作用亜群へ制限した `Φ`。 -/
abbrev AgPhiAct : MonoidOn.{0, 0, 0} DAct :=
  AgPhi.restrict (ActionCategory.π Gam ℕ) dact_fsm

/-- ★★底圏の関手に沿って `𝔽` を押し出す関手。

★`Div` と `deg` はそのまま、`base` だけ `G` で送る。
★**`Φ.restrict` の `map` は `Φ.map (G.map ·)` そのものなので、`div` 成分は `rfl`。** -/
def elemPush {D : Type u} [Category.{v} D] {D' : Type u3} [Category.{v3} D']
    (Φ : MonoidOn.{v, u, w} D) (G : D' ⥤ D)
    (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α)) :
    ElemFrobCat (Φ.restrict G hG) ⥤ ElemFrobCat Φ where
  obj W := ⟨G.obj W.base⟩
  map {X Y} f := { base := G.map f.base, div := f.div, deg := f.deg }
  map_id W := ElemFrobCat.Hom.ext (by simp) rfl rfl
  map_comp f g := ElemFrobCat.Hom.ext (by simp) rfl rfl

/-- ★反例の圏 —— 作用亜群上の `𝔽`。 -/
abbrev AgC : Type := ElemFrobCat AgPhiAct

/-- ★★★**`𝒟 = SingleObj Γ` 上の pre-Frobenioid**。

★★**ここが要点** —— 圏は作用亜群の上の `𝔽` だが、
★**`Base` は 1 対象の底圏に落ちる。** -/
def AgP : PreFrobenioid AgC AgPhi where
  toElem := elemPush AgPhi (ActionCategory.π Gam ℕ) dact_fsm
  divisorial _ := isDivisorial_nat
  totEpiC := isTotallyEpimorphic_elemFrobCat dact_totEpi (fun _ => isIntegralMonoid_nat)
  totEpiD := dsingle_totEpi
  connectedC := isConnected_elemFrobCat AgPhiAct
  connectedD := dsingle_isConnected

/-- ★`DAct` の射の**台となる群の元**。 -/
def dactElt {X Y : DAct} (f : X ⟶ Y) : Gam := f.1

theorem dactElt_spec {X Y : DAct} (f : X ⟶ Y) : dactElt f • X.back = Y.back := f.2

/-! ## ★★2 つの事実 —— base-trivial なのに Aut-ample でない -/

/-- ★★**`𝒞` の対象はすべて同型** —— 作用が推移的だから。 -/
theorem ag_allIso (X Y : AgC) : Nonempty (X ≅ Y) := by
  obtain ⟨γ, hγ⟩ := gam_transitive X.base.back Y.base.back
  let b : X.base ⟶ Y.base := ⟨γ, hγ⟩
  haveI hbi : IsIso b := IsIso.of_groupoid b
  let f : X ⟶ Y := { base := b, div := (0 : ℕ), deg := 1 }
  haveI : IsIso f := (ElemFrobCat.isIso_iff f).mpr ⟨hbi, isAddUnit_zero, rfl⟩
  exact ⟨asIso f⟩

/-- ★★★**すべての対象が base-trivial** ——
`𝒟` が 1 対象なので `BaseIsomorphic` は常に真で、対象がすべて同型だから。 -/
theorem ag_baseTrivial (X : AgC) : IsBaseTrivial AgP X :=
  fun Dd _ => ag_allIso Dd X

/-- ★反例の対象 —— `0 ∈ ℕ` の上。 -/
abbrev AgA : AgC := ⟨(⟨(), (0 : ℕ)⟩ : DAct)⟩

/-- ★★★**それは `Aut-ample` ではない** ——
`Aut_𝒞(A)` の像は `Stab(0)` であって `Γ` ではない。

★★**したがって `Proposition 1.6, (v)` の `⟸` の証明に要る
「底を指定した同型の取り直し」は、この `𝒞` では効かない。** -/
theorem ag_not_autAmple : ¬ IsAutAmple AgP AgA := by
  intro h
  obtain ⟨γ₀, hγ₀⟩ := gam_stab_proper
  obtain ⟨φ, -, hbase⟩ := h γ₀ (IsIso.of_groupoid _)
  -- ★`φ.base` は `DAct` の射なので、その台元は `0` を固定する
  have hfix := dactElt_spec (ElemFrobCat.Hom.base φ)
  have hb : dactElt (ElemFrobCat.Hom.base φ) = γ₀ := hbase
  rw [hb] at hfix
  exact hγ₀ hfix

/-! ## ★★★設計は `Definition 1.3` を満たさない —— 自己訂正

★**上の設計は `preStepSpan`(`Definition 1.3, (i), (b)`)を満たさない。**
★★**それをここで式にする** —— 散文の訂正ではなく、証明で示す。 -/

/-- ★★★**この `𝒞` は `preStepSpan` を満たさない**。

★★**要点**: pre-step `φ : X ⟶ A` の底は「`X` の点を `A` の点へ送る」`Γ` の元に
限られる。したがって `inv(Base φ) ≫ Base ψ` は必ず `A` の点を `B` の点へ送り、
★**`Γ` 全体を張れない。**

★**したがって `Gap_1_6_v` の反例としては、この設計は使えない。**
★★**`preStepSpan` は「底の同型を pre-step の比で実現せよ」と要求しており、
それ自体が `Aut-ample` に近い条件である** ——
これは私が最初の見立てで**使っていなかった道具**である。 -/
theorem ag_not_frobenioidCore : ¬ FrobenioidCore AgP := by
  intro F
  obtain ⟨γ₀, hγ₀⟩ := gam_stab_proper
  obtain ⟨X, φ, ψ, hφ, -, heq⟩ := F.preStepSpan AgA AgA γ₀ (IsIso.of_groupoid _)
  have hφs : dactElt (ElemFrobCat.Hom.base φ) • X.base.back = (0 : ℕ) :=
    dactElt_spec (ElemFrobCat.Hom.base φ)
  have hψs : dactElt (ElemFrobCat.Hom.base ψ) • X.base.back = (0 : ℕ) :=
    dactElt_spec (ElemFrobCat.Hom.base ψ)
  -- ★`inv` を消す —— 両辺に `Base φ` を前から掛ければ `hom_inv_id` で潰れる
  have h2 := congrArg (fun t : (AgP.toElem.obj AgA).base ⟶ (AgP.toElem.obj AgA).base =>
    AgP.Base φ ≫ t) heq
  rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp] at h2
  -- ★`SingleObj` では `f ≫ g = g * f` なので、これは群の等式そのもの
  have hkey : γ₀ * dactElt (ElemFrobCat.Hom.base φ) = dactElt (ElemFrobCat.Hom.base ψ) := h2
  refine hγ₀ ?_
  have h3 : (γ₀ * dactElt (ElemFrobCat.Hom.base φ)) • X.base.back = (0 : ℕ) := by
    rw [hkey]; exact hψs
  rw [mul_smul, hφs] at h3
  exact h3

end ABC3.Check.FrdI
