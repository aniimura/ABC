import ABC3.Found.PGC.StageUpperNaturality
import ABC3.Found.PGC.QpResidueField

/-!
# [pGC] Theorem 4.2 の全単射性 —— 単射性を「純群論の 1 文」と**同値**にする

`Skeleton/PGC/Section4.lean::theorem_4_2` は 2026-09-08 に仮説が 0 本になり

```
theorem theorem_4_2 (K K' : PAdicLocalField p) :
    Function.Bijective
      (naturalOuterIso (ABC3.Found.PGC.ramificationFiltration p)
        (ABC3.Found.PGC.isNaturalFiltration_ramificationFiltration p) (K := K) (K' := K'))
```

という形になった。本ファイルはその **全単射性の両側**を扱う。

原文 (pGC 物理 p.7, Theorem 4.2 の証明):

> Proof: That this morphism is injective follows by looking at the induced isomorphism
> between Γab K ∼= (K×)∧and Γab K′ ∼= ((K′)×)∧. (Alternatively, the more group-theoretically
> oriented reader may prefer to regard the injectivity of this morphism as a consequence of
> the fact that the centralizer of ΓK in ΓQp is trivial.) Now let us show that the morphism
> is surjective. Let α : ΓK ∼= ΓK′ be an isomorphism that is compatible with the filtrations.
> Then, in the preceding paragraph, we constructed a field isomorphism αK : K ∼= K′. To
> see that the morphism induced by αK is equal to the original α (up to composition with
> an inner isomorphism), it suffices simply to construct the analogues αL : L ∼= L′ of αK for
> corresponding finite extensions L and L′ of K and K′, and then conclude via a standard
> general nonsense argument.

## 何を得たか（測ってあること・測っていないこと）

| 段 | 結果 |
|---|---|
| 単射性 | ★`CentralizerActsTriviallyOnBase K` と **同値**（`injective_naturalOuterIso_iff`）。 |
| 全射性 | ★「濾過つき同型が代数閉包の環同型による共役で書ける」に**還元**した（`surjective_naturalOuterIso_of_forall_extension`）。★中身（Cor 3.3・Lemma 4.1・general nonsense）は**未着手**。 |
| 全単射 | ★上の 2 つを合わせた `bijective_naturalOuterIso_ramificationFiltration`。 |

★**単射性は「無条件に閉じた」のではない。**★`CentralizerActsTriviallyOnBase` は
**未証明の 1 文**であり、本ファイルはそれが単射性と**過不足なく同じ**であることを示しただけである。
★ただしこれは原典の括弧書きが主張していることそのものであり、
★**同値**まで言えたので「強すぎる仮定を置いた」ということは起こり得ない。

## ★原典との差（逸脱の記録 1）

原典は「the centralizer of Γ_K in Γ_Qp is trivial」と書く。本ファイルが使うのは
それより**弱い**

  C_{Γ_Qp}(Γ_K) ⊆ Γ_K   （＝「Γ_K を中心化する元は K の上で恒等に働く」）

である（`CentralizerActsTriviallyOnBase`）。★弱い方で十分であるうえ、
★**弱い方は単射性と同値**（`injective_naturalOuterIso_iff`）——すなわち
原典の「trivial」は**必要以上に強い**。これが本ファイルで見つかった短い道である。

## ★原典の 2 本のうち、こちらを採った理由（測定つき）

原典が先に挙げる道は `Γ_K^ab ≅ (K^×)^` を経由するもので、これには
**Artin 写像が体の同型と同変であること**が要る。本リポジトリではそれは
`Found/PGC/ArtinEquivariance.lean::ArtinEquivariance`（および同値な 4 形）として
**まだ開いている壁**である。★括弧書きの群論的な道は、その壁を 1 度も踏まない。

## ★抽象核（分岐・付値・Galois の語彙が 1 語も出ない）

`exists_comm_of_conj_eq_conj` は
「2 つの環同型 `s t : A ≃+* B` による共役が、`B` の内部同型のぶんだけ違うなら、
`t = c ∘ s ∘ τ` と書けて `τ` は問題の族と可換」
という**一般の半環についての 1 文**である。分岐も付値も Galois も出てこない。
具体層はここに `A := K̄`, `B := K̄′`, 族 `:= Γ_K` を代入するだけ。

## ★退化していないことの確認（`Check/PGC/Theorem42NaiveGC.lean` への回答）

`Check/PGC/Theorem42NaiveGC.lean` は「`∀ RF` の形だと `RF := topFiltration p` を
代入できて素朴 Grothendieck 予想が出てしまう」と示した。現行形は `RF` を実物
`ramificationFiltration p` に固定しており、
`Found/PGC/Prop22FilteredHypothesis.lean::ramificationFiltration_Gv_zero_ne_top`
（無条件）が `Γ_K^0 ≠ ⊤` を与えるので、実物は `topFiltration` **ではない**。
また `::map_ramificationFiltration_iff_pos` が「全 `v` で保つ」と
「`v > 0` で保つ」が同値だと言うので、原典が `index > 0` に限っていることとも
食い違わない。★すなわち現行形は退化しても強すぎもしていない。

## ★測っていないこと（正直に）

* `Γ_K^v`（`v > 0`）が `v` について実際に真に減っていくかは**測っていない**。
  上の 2 点は「`v = 0` で `⊤` でない」「`v ≤ 0` は無内容」までしか言わない。
* 全射性の中身（Cor 3.3 → `α_K` → general nonsense）は**着手していない**。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

/-! ## §1 抽象核 —— 一般の半環の同型についての 1 文

★ここには分岐・付値・Galois の語彙が 1 語も無い。 -/

section AbstractCore

variable {A B : Type*} [Semiring A] [Semiring B]

/-- ★★★**抽象核**。`s t : A ≃+* B` を 2 つの同型、`c : B ≃+* B` とする。
`A` の自己同型の族 `G : ι → (A ≃+* A)` について

  `t ∘ G i ∘ t⁻¹ = c ∘ (s ∘ G i ∘ s⁻¹) ∘ c⁻¹`   （すべての `i`）

が成り立つなら、`τ := s⁻¹ ∘ c⁻¹ ∘ t` は

* `c ∘ s ∘ τ = t`（`t` は `c ∘ s` と `τ` のぶんだけしか違わない）、
* すべての `G i` と**可換**

を満たす。★「共役が内部同型のぶんだけ一致する ⟹ もとの同型は中心化元のぶんだけ違う」
という、torsor の言葉を使わずに書いた形である。証明は 3 行。 -/
theorem exists_comm_of_conj_eq_conj {ι : Type*} (s t : A ≃+* B) (c : B ≃+* B)
    (G : ι → (A ≃+* A))
    (h : ∀ (i : ι) (x : B), t (G i (t.symm x)) = c (s (G i (s.symm (c.symm x))))) :
    ∃ τ : A ≃+* A,
      (∀ y : A, c (s (τ y)) = t y) ∧ (∀ (i : ι) (y : A), τ (G i y) = G i (τ y)) := by
  refine ⟨t.trans (c.symm.trans s.symm), ?_, ?_⟩
  · intro y
    show c (s (s.symm (c.symm (t y)))) = t y
    simp
  · intro i y
    show s.symm (c.symm (t (G i y))) = G i (s.symm (c.symm (t y)))
    have hx := h i (t y)
    rw [RingEquiv.symm_apply_apply] at hx
    rw [hx, RingEquiv.symm_apply_apply, RingEquiv.symm_apply_apply]

end AbstractCore

variable {p : ℕ} [Fact p.Prime]

/-! ## §2 延長の選択非依存性の**一様**版

`GaloisTransfer.lean::galMulEquivOf_indep` は `∃ c` が `∀ g` の**内側**にある。
`FilteredGroup.Iso.setoid` は `g` によらない 1 つの `c` を要求するので、
そのままでは使えない。★証明本体は `c` を `g` に依らず作っているので、
束縛の順を入れ替えるだけで一様版になる。 -/

/-- `galMulEquivOf_indep` の一様版（`∃ c` を `∀ g` の外に出したもの）。

★証明は `Found/PGC/GaloisTransfer.lean::galMulEquivOf_indep` と同一で、
`refine` の位置だけが違う。 -/
theorem galMulEquivOf_indep_uniform {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ1 ᾱ2 : K.closure ≃+* K'.closure)
    (hfwd1 : ∀ x : K.carrier, ᾱ1 (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (hfwd2 : ∀ x : K.carrier, ᾱ2 (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    ∃ c : K'.absGal, ∀ g : K.absGal,
      galMulEquivOf α ᾱ1 hfwd1 g = c * galMulEquivOf α ᾱ2 hfwd2 g * c⁻¹ := by
  set c0 : K.closure ≃+* K.closure := ᾱ1.trans ᾱ2.symm with hc0
  set c : K.absGal := AlgEquiv.ofRingEquiv (f := c0) (galMulEquivOf_indep_fix α ᾱ1 ᾱ2 hfwd1 hfwd2)
    with hc
  have hc_apply : ∀ z : K.closure, c z = ᾱ2.symm (ᾱ1 z) := fun z => rfl
  have hcinv_apply : ∀ z : K.closure, c⁻¹ z = ᾱ1.symm (ᾱ2 z) := by
    intro z
    have hcc : c (c⁻¹ z) = z := by
      show (c * c⁻¹) z = z
      rw [mul_inv_cancel]; rfl
    rw [hc_apply] at hcc
    have h2 : ᾱ1 (c⁻¹ z) = ᾱ2 z := by
      have := congrArg ᾱ2 hcc
      rwa [RingEquiv.apply_symm_apply] at this
    have := congrArg ᾱ1.symm h2
    rwa [RingEquiv.symm_apply_apply] at this
  refine ⟨galMulEquivOf α ᾱ2 hfwd2 c, fun g => ?_⟩
  apply AlgEquiv.ext
  intro x
  show ᾱ1 (g.toRingEquiv (ᾱ1.symm x))
      = (galMulEquivOf α ᾱ2 hfwd2 c * galMulEquivOf α ᾱ2 hfwd2 g
          * (galMulEquivOf α ᾱ2 hfwd2 c)⁻¹) x
  rw [← map_inv, ← map_mul, ← map_mul]
  show ᾱ1 (g.toRingEquiv (ᾱ1.symm x)) = ᾱ2 ((c * g * c⁻¹).toRingEquiv (ᾱ2.symm x))
  have hmul : (c * g * c⁻¹).toRingEquiv (ᾱ2.symm x) = c (g (c⁻¹ (ᾱ2.symm x))) := rfl
  rw [hmul, hcinv_apply, RingEquiv.apply_symm_apply, hc_apply, RingEquiv.apply_symm_apply]
  rfl

/-- 既定の延長 `extendToClosure α` を基準にした一様版。 -/
theorem galMulEquiv_conj_indep_uniform {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) :
    ∃ c : K'.absGal, ∀ g : K.absGal,
      galMulEquivOf α ᾱ hfwd g = c * galMulEquiv α g * c⁻¹ :=
  galMulEquivOf_indep_uniform α ᾱ (extendToClosure α) hfwd (extendToClosure_algebraMap α)

/-! ## §3 恒等の濾過つき同型 -/

/-- 恒等写像は濾過つき同型（`Subgroup.map id = id` だけ）。 -/
def reflFilteredIso (A : FilteredGroup) : FilteredGroup.Iso A A where
  equiv := ContinuousMulEquiv.refl A.G
  map_Gv v := by
    ext g
    constructor
    · rintro ⟨h, hh, rfl⟩; exact hh
    · intro hg; exact ⟨g, hg, rfl⟩

/-! ## §4 外部同型類の同定 —— 全射性はここに還元される -/

/-- 濾過つき同型 `F` が「`K̄` から `K̄′` への環同型 `ᾱ` による共役」で書けているなら、
`F` の外部同型類は `naturalOuterIso RF hnat α` に**等しい**。

★`ᾱ` は `extendToClosure α` である必要はない（`galMulEquiv_conj_indep_uniform`
が延長の差を吸収する）。★また「内部同型のぶんだけ違う `ᾱ`」を許しても一般性は
増えない —— `c ∘ ᾱ`（`c ∈ Γ_{K′}`）もまた `α` と互換な延長だからである。 -/
theorem naturalOuterIso_eq_mk_of_conj (RF : RamificationFiltration p)
    (hnat : IsNaturalFiltration RF) {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (F : FilteredGroup.Iso (filtOf RF K) (filtOf RF K'))
    (hF : ∀ g : K.absGal, (F.equiv g : K'.absGal) = galMulEquivOf α ᾱ hfwd g) :
    naturalOuterIso RF hnat α = Quotient.mk _ F := by
  obtain ⟨c, hc⟩ := galMulEquiv_conj_indep_uniform α ᾱ hfwd
  refine Quotient.sound ⟨c, ?_⟩
  intro g
  rw [hF g]
  exact hc g

def naturalOuterIso_eq_mk_of_conj.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-- `hF` を点ごとの形で与えるための橋渡し。 -/
theorem eq_galMulEquivOf_of_apply {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x))
    (u : K'.absGal) (g : K.absGal)
    (h : ∀ x : K'.closure, u x = ᾱ (g (ᾱ.symm x))) :
    u = galMulEquivOf α ᾱ hfwd g :=
  AlgEquiv.ext h

/-- ★★★**全射性の還元**。

「どの濾過つき同型 `F : Γ_K ≅ Γ_{K′}` も、`K` を `K′` に写す `K̄ ≃+* K̄′` による
共役で書ける」が示せれば、`naturalOuterIso` は全射である。

★これが原典の証明の残り全部（Corollary 3.3 で `α_K : K ≅ K′` を作り、
Lemma 4.1 で `I` が `K` を張ることを使い、有限次拡大 `L` について `α_L` を作って
"standard general nonsense argument" で貼り合わせる）に対応する。
★**本ファイルはその中身には手を付けていない。** -/
theorem surjective_naturalOuterIso_of_forall_extension (RF : RamificationFiltration p)
    (hnat : IsNaturalFiltration RF) {K K' : PAdicLocalField p}
    (h : ∀ F : FilteredGroup.Iso (filtOf RF K) (filtOf RF K'),
      ∃ (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
        (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
          = algebraMap K'.carrier K'.closure (α x)),
        ∀ g : K.absGal, (F.equiv g : K'.absGal) = galMulEquivOf α ᾱ hfwd g) :
    Function.Surjective (naturalOuterIso RF hnat (K := K) (K' := K')) := by
  intro Φ
  obtain ⟨F, rfl⟩ := Quotient.exists_rep Φ
  obtain ⟨α, ᾱ, hfwd, hF⟩ := h F
  exact ⟨α, naturalOuterIso_eq_mk_of_conj RF hnat α ᾱ hfwd F hF⟩

def surjective_naturalOuterIso_of_forall_extension.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 8, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-! ## §5 中心化群の条件 —— 単射性と**同値** -/

/-- ★★原典 p.7 の括弧書き「the centralizer of Γ_K in Γ_Qp is trivial」の**弱形**。

`Γ_Qp = Aut(K̄/Q_p)` の元 `σ`（`K` の上で `ρ` として働くもの）が `Γ_K` を
中心化するなら `ρ = id`、すなわち `C_{Γ_Qp}(Γ_K) ⊆ Γ_K`。

★原典の「trivial」（`C = 1`）より弱い。★それでも単射性には十分であり、
しかも**同値**である（`injective_naturalOuterIso_iff`）。

★`σ` が `ℚ_p` を固定することは仮定していない（`ρ` が `ℚ_p`-代数同型であることから
自動的に従うので、仮定に入れると主張が弱くなるだけである）。 -/
def CentralizerActsTriviallyOnBase (K : PAdicLocalField p) : Prop :=
  ∀ (ρ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (σ : K.closure ≃+* K.closure),
    (∀ x : K.carrier, σ (algebraMap K.carrier K.closure x)
      = algebraMap K.carrier K.closure (ρ x)) →
    (∀ (g : K.absGal) (x : K.closure), σ (g x) = g (σ x)) →
    ρ = AlgEquiv.refl

def CentralizerActsTriviallyOnBase.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-- ★★★★**Theorem 4.2 の単射性**（中心化群の条件から）。

`naturalOuterIso RF hnat α = naturalOuterIso RF hnat β` は
「`β̄` による共役と `ᾱ` による共役が `Γ_{K′}` の内部同型のぶんだけ一致する」ことなので、
抽象核 `exists_comm_of_conj_eq_conj` が `τ = ᾱ⁻¹ ∘ c⁻¹ ∘ β̄` を返す。
`τ` は `Γ_K` を中心化し、`K` の上では `α⁻¹ ∘ β` として働くので、
仮説から `α⁻¹ ∘ β = id`、すなわち `α = β`。

★★**濾過を 1 度も使っていない** —— `RF` は何でもよい。
原典が `Out_Filt` で述べていても、単射性の側に濾過は要らない。 -/
theorem injective_naturalOuterIso_of_centralizer (RF : RamificationFiltration p)
    (hnat : IsNaturalFiltration RF) {K K' : PAdicLocalField p}
    (hC : CentralizerActsTriviallyOnBase K) :
    Function.Injective (naturalOuterIso RF hnat (K := K) (K' := K')) := by
  intro α β hab
  obtain ⟨c₀, hc⟩ := Quotient.exact hab
  set c : K'.absGal := c₀ with hcdef
  have hconj : ∀ (g : K.absGal) (x : K'.closure),
      extendToClosure β (g ((extendToClosure β).symm x))
        = c (extendToClosure α (g ((extendToClosure α).symm (c.symm x)))) :=
    fun g x => congrArg (fun u : K'.absGal => u x) (hc g)
  obtain ⟨τ, hτ1, hτ2⟩ :=
    exists_comm_of_conj_eq_conj (extendToClosure α) (extendToClosure β) c.toRingEquiv
      (fun g : K.absGal => g.toRingEquiv) hconj
  have hcomm : ∀ (g : K.absGal) (y : K.closure), τ (g y) = g (τ y) := hτ2
  have hfwdτ : ∀ x : K.carrier,
      τ (algebraMap K.carrier K.closure x)
        = algebraMap K.carrier K.closure ((β.trans α.symm) x) := by
    intro x
    have hL : c.toRingEquiv (extendToClosure α (τ (algebraMap K.carrier K.closure x)))
        = algebraMap K'.carrier K'.closure (β x) := by
      rw [hτ1 (algebraMap K.carrier K.closure x)]
      exact extendToClosure_algebraMap β x
    have hR : c.toRingEquiv (extendToClosure α (algebraMap K.carrier K.closure (α.symm (β x))))
        = algebraMap K'.carrier K'.closure (β x) := by
      show c (extendToClosure α (algebraMap K.carrier K.closure (α.symm (β x))))
        = algebraMap K'.carrier K'.closure (β x)
      rw [extendToClosure_algebraMap α, AlgEquiv.apply_symm_apply]
      exact c.commutes (β x)
    exact (extendToClosure α).injective (c.toRingEquiv.injective (hL.trans hR.symm))
  have hρ := hC (β.trans α.symm) τ hfwdτ hcomm
  refine AlgEquiv.ext fun x => ?_
  have h0 : α.symm (β x) = x :=
    congrArg (fun u : K.carrier ≃ₐ[ℚ_[p]] K.carrier => u x) hρ
  have h1 := congrArg α h0
  rw [AlgEquiv.apply_symm_apply] at h1
  exact h1.symm

def injective_naturalOuterIso_of_centralizer.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-- ★★**逆向き** —— 単射性から中心化群の条件が出る。

`σ` が `Γ_K` を中心化するなら、`σ` を延長に取った共役は `Γ_K` 上の**恒等**である。
ゆえに `ρ` の外部同型類と `id` の外部同型類が一致し、単射性から `ρ = id`。

★これで仮説 `CentralizerActsTriviallyOnBase` は「強すぎる」ことがあり得ない。 -/
theorem centralizerActsTriviallyOnBase_of_injective (RF : RamificationFiltration p)
    (hnat : IsNaturalFiltration RF) {K : PAdicLocalField p}
    (hinj : Function.Injective (naturalOuterIso RF hnat (K := K) (K' := K))) :
    CentralizerActsTriviallyOnBase K := by
  intro ρ σ hfwd hcomm
  have h1 : naturalOuterIso RF hnat ρ = Quotient.mk _ (reflFilteredIso (filtOf RF K)) :=
    naturalOuterIso_eq_mk_of_conj RF hnat ρ σ hfwd _ (fun g => AlgEquiv.ext fun x => by
      show g x = σ (g (σ.symm x))
      rw [hcomm g (σ.symm x), RingEquiv.apply_symm_apply])
  have h2 : naturalOuterIso RF hnat (AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier)
      = Quotient.mk _ (reflFilteredIso (filtOf RF K)) :=
    naturalOuterIso_eq_mk_of_conj RF hnat (AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier)
      (RingEquiv.refl K.closure) (fun _ => rfl) _ (fun _ => AlgEquiv.ext fun _ => rfl)
  exact hinj (h1.trans h2.symm)

/-- ★★★★★**単射性 ⟺ 中心化群の条件**。

★原典が括弧で「こう思ってもよい」と書いた言い換えが、実は**同値**であることを言う。
★左辺は `K′` を `K` に取った場合だが、右辺からは**すべての `K′`** について
単射性が出る（`injective_naturalOuterIso_of_centralizer`）。 -/
theorem injective_naturalOuterIso_iff (RF : RamificationFiltration p)
    (hnat : IsNaturalFiltration RF) (K : PAdicLocalField p) :
    Function.Injective (naturalOuterIso RF hnat (K := K) (K' := K))
      ↔ CentralizerActsTriviallyOnBase K :=
  ⟨centralizerActsTriviallyOnBase_of_injective RF hnat,
   fun hC => injective_naturalOuterIso_of_centralizer RF hnat hC⟩

def injective_naturalOuterIso_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-! ## §6 非空虚性（★弱い witness であることを明記する） -/

/-- `Aut(K/ℚ_p)` が自明なら条件は自明に成り立つ。 -/
theorem centralizerActsTriviallyOnBase_of_subsingleton (K : PAdicLocalField p)
    [Subsingleton (K.carrier ≃ₐ[ℚ_[p]] K.carrier)] :
    CentralizerActsTriviallyOnBase K :=
  fun _ _ _ _ => Subsingleton.elim _ _

/-- `ℚ_p` の `ℚ_p`-代数自己同型は恒等だけ。 -/
theorem subsingleton_selfAlgEquiv (p : ℕ) [Fact p.Prime] :
    Subsingleton (ℚ_[p] ≃ₐ[ℚ_[p]] ℚ_[p]) :=
  ⟨fun e f => AlgEquiv.ext fun x => by
    have he : e x = x := by simpa using e.commutes x
    have hf : f x = x := by simpa using f.commutes x
    rw [he, hf]⟩

/-- ★**非空虚性の witness**: `K = ℚ_p` では条件が成り立つ。

☆★**弱い witness であることを隠さない**: `Aut(ℚ_p/ℚ_p)` が 1 点なので、
この `K` については単射性そのものが定義域の subsingleton 性から自明に従う。
★したがってこの witness は「条件が矛盾していない」ことしか言わない。
★**条件が過不足ないこと**の根拠は `injective_naturalOuterIso_iff` の方である。 -/
theorem centralizerActsTriviallyOnBase_selfField (p : ℕ) [Fact p.Prime] :
    CentralizerActsTriviallyOnBase (selfField p) := by
  haveI : Subsingleton ((selfField p).carrier ≃ₐ[ℚ_[p]] (selfField p).carrier) :=
    subsingleton_selfAlgEquiv p
  exact centralizerActsTriviallyOnBase_of_subsingleton _

/-! ## §7 実物の濾過（`theorem_4_2` の形）への特殊化 -/

/-- ★実物の高次分岐濾過での単射性。 -/
theorem injective_naturalOuterIso_ramificationFiltration {K K' : PAdicLocalField p}
    (hC : CentralizerActsTriviallyOnBase K) :
    Function.Injective (naturalOuterIso (ramificationFiltration p)
      (isNaturalFiltration_ramificationFiltration p) (K := K) (K' := K')) :=
  injective_naturalOuterIso_of_centralizer _ _ hC

/-- `naturalOuterIsoOfAlgEquiv`（`StageUpperNaturality.lean` の仮定なし版）での単射性。 -/
theorem injective_naturalOuterIsoOfAlgEquiv {K K' : PAdicLocalField p}
    (hC : CentralizerActsTriviallyOnBase K) :
    Function.Injective (naturalOuterIsoOfAlgEquiv (p := p) (K := K) (K' := K')) :=
  injective_naturalOuterIso_ramificationFiltration hC

/-- ★★★★★★**[pGC] Theorem 4.2 の全単射性が、残り 2 つの文に還元された**。

★`Skeleton/PGC/Section4.lean::theorem_4_2` の結論とちょうど同じ形。
★残っているのは
(1) `CentralizerActsTriviallyOnBase K`（★単射性と**同値**なので、これ以上弱められない）、
(2) 濾過つき同型が代数閉包の環同型による共役で書けること(★原典の Cor 3.3 + Lemma 4.1 +
    "standard general nonsense argument" に対応)
の 2 つだけである。 -/
theorem bijective_naturalOuterIso_ramificationFiltration {K K' : PAdicLocalField p}
    (hC : CentralizerActsTriviallyOnBase K)
    (hS : ∀ F : FilteredGroup.Iso (filtOf (ramificationFiltration p) K)
        (filtOf (ramificationFiltration p) K'),
      ∃ (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
        (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
          = algebraMap K'.carrier K'.closure (α x)),
        ∀ g : K.absGal, (F.equiv g : K'.absGal) = galMulEquivOf α ᾱ hfwd g) :
    Function.Bijective (naturalOuterIso (ramificationFiltration p)
      (isNaturalFiltration_ramificationFiltration p) (K := K) (K' := K')) :=
  ⟨injective_naturalOuterIso_ramificationFiltration hC,
   surjective_naturalOuterIso_of_forall_extension _ _ hS⟩

def bijective_naturalOuterIso_ramificationFiltration.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

#print axioms exists_comm_of_conj_eq_conj
#print axioms galMulEquivOf_indep_uniform
#print axioms galMulEquiv_conj_indep_uniform
#print axioms reflFilteredIso
#print axioms naturalOuterIso_eq_mk_of_conj
#print axioms eq_galMulEquivOf_of_apply
#print axioms surjective_naturalOuterIso_of_forall_extension
#print axioms injective_naturalOuterIso_of_centralizer
#print axioms centralizerActsTriviallyOnBase_of_injective
#print axioms injective_naturalOuterIso_iff
#print axioms centralizerActsTriviallyOnBase_of_subsingleton
#print axioms subsingleton_selfAlgEquiv
#print axioms centralizerActsTriviallyOnBase_selfField
#print axioms injective_naturalOuterIso_ramificationFiltration
#print axioms injective_naturalOuterIsoOfAlgEquiv
#print axioms bijective_naturalOuterIso_ramificationFiltration

end ABC3.Found.PGC
