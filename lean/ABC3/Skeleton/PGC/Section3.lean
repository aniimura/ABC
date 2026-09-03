import ABC3.Skeleton.PGC.Section2

/-!
# [pGC] §3 — 命題

設定・記号は `ABC3/Skeleton/PGC/Setup.lean`。§2 の `RamificationFiltration`・
`RecoverableAsAddModule` を継続して使う。

構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic
Local Fields/section-3.html`(PDF 目視確認 2026-09-03、物理 p.6)。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## Corollary 3.1 -/

/-- **[pGC] Corollary 3.1**

原文 (pGC p.6):
> Given a continuous Q_p[Γ_K]-vector space of finite Q_p-dimension, the issue of whether or
> not V is Hodge-Tate (as well as the invariants d_V(i)) can be determined entirely
> group-theoretically from the filtered group Γ_K.

## 形式化上の簡略化(逸脱として記録)

「Hodge-Tate である」の定義(`d_V ≤ dim_{Q_p}(V)`、`d_V(i)` は `V(−i)⊗K̄^` の Γ_K-不変部分の
次元)は Proposition 2.2 の `CompKbar`(K̄^ の抽象化)とテンソル積・Tate twist を要し、
本セッションでは未構築。ここでは「Hodge-Tate である」を**抽象的な述語パラメータ**
`isHodgeTate : ∀ K, Prop` として受け取り、その値が α で保たれることだけを主張する
——実際の Hodge-Tate 述語の構成は独立した課題として残す。

## 依拠する境界外の結果

- [1] Serre, *Abelian ℓ-adic Representations*, Chapter III §1.2(d_V ≤ dim_{Q_p}(V)、
  Hodge-Tate の定義そのもの)。mathlib 不在。
- Proposition 2.2(K̄^ の回復)への直接依存。 -/
theorem cor_3_1 (_RF : RamificationFiltration p)
    (V : PAdicLocalField p → Type*) [∀ K, AddCommGroup (V K)] [∀ K, Module ℚ_[p] (V K)]
    [∀ K, FiniteDimensional ℚ_[p] (V K)] [∀ K, SMul K.absGal (V K)]
    (isHodgeTate : ∀ K, Prop) :
    ∀ {K K' : PAdicLocalField p} (_α : ContinuousMulEquiv K.absGal K'.absGal),
      isHodgeTate K ↔ isHodgeTate K' := sorry

def cor_3_1.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def cor_3_1.needs : List ProofObligation :=
  [ .citation "[1] Serre, Abelian ℓ-adic Representations, Chapter III §1.2"
      "Hodge-Tate の定義(d_V ≤ dim_{Q_p}(V))"
      (.absent "mathlib v4.31.0-rc2 実測: Hodge-Tate 表現論に相当する宣言はゼロ件") 6,
    .otherPaper "pGC" "Proposition 2.2" 5 ]

/-! ## Definition 3.2 -/

/-- **[pGC] Definition 3.2** — uniformizing な E[Γ_K]-加群。`sorry` 無し(純粋な定義)。

原文 (pGC p.6):
> We shall call the E[Γ_K]-module V uniformizing if the restriction of ρ_V to some open
> subgroup I of U_K (⊆ Γ^a_K^b) is the morphism I → E× induced by restricting some morphism
> of fields K → E to I ⊆ U_K ⊆ K.

## 形式化上の簡略化(逸脱として記録)

原文の「ρ_V は U_K ⊆ Γ_K^ab 上に制限できる」は、**局所類体論の相互律**
Γ_K^ab ≅ (K^×)^(§1 Proposition 1.2 の論拠、mathlib 不在)を暗黙に経由する——
U_K は本来 K^× の部分群であり、Γ_K の部分群ではない。

ここでは相互律を明示的なパラメータ `toGal : U_K → Γ_K`(未構築の辞書、まだ本物には
なっていない)として受け取ることで、定義全体を条件付きで well-typed にする
——`toGal` が本物になれば(Interface が実装されれば)この定義はそのまま使える。 -/
def IsUniformizing (K : PAdicLocalField p) (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) (ρ : K.absGal →* Eˣ) : Prop :=
  ∃ (I : Set K.carrier) (hIU : I ⊆ {x : K.carrier | ‖x‖ = 1}) (_hopen : IsOpen I)
    (ι : K.carrier →+* E), ∀ x (hx : x ∈ I), ((ρ (toGal ⟨x, hIU hx⟩) : Eˣ) : E) = ι x

def IsUniformizing.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Definition 3.2", sectionId := "def-3-2" }

/-! ## Corollary 3.3 -/

/-- **[pGC] Corollary 3.3**

原文 (pGC p.6):
> Given a continuous E[Γ_K]-module V of E-dimension 1, the issue of whether or not V is
> uniformizing can be determined entirely group-theoretically from the filtered group Γ_K.

## 依拠する境界外の結果

- [1] Serre, Chapter III, Appendix §5(d_V(i) による uniformizing の判定条件式)。
- Corollary 3.1 への直接依存(d_V(i) が回復できることを使う)。 -/
theorem cor_3_3 (_RF : RamificationFiltration p)
    (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : ∀ K : PAdicLocalField p, {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) :
    ∀ {K K' : PAdicLocalField p} (_α : ContinuousMulEquiv K.absGal K'.absGal)
      (ρ : K.absGal →* Eˣ) (ρ' : K'.absGal →* Eˣ),
      IsUniformizing K E (toGal K) ρ ↔ IsUniformizing K' E (toGal K') ρ' := sorry

def cor_3_3.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.3", sectionId := "cor-3-3" }

def cor_3_3.needs : List ProofObligation :=
  [ .citation "[1] Serre, Abelian ℓ-adic Representations, Chapter III, Appendix §5"
      "d_V(1) = [E:K]; d_V(0) = [E:K]·([K:Q_p]−1) による uniformizing の判定"
      (.absent "mathlib v4.31.0-rc2 実測: 該当なし") 6,
    .otherPaper "pGC" "Corollary 3.1" 6 ]

end ABC3.Skeleton.PGC
