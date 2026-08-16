import ABC3.Found.FrdI.WitnessFrobenioid

/-!
# ★負の対照 —— **Frobenioid でない pre-Frobenioid**

`WitnessFrobenioid.lean` は「**`Frobenioid` は空でない**」を示した
(`wIsFrobenioid : Frobenioid wP`)。しかしそれだけでは

> `Definition 1.3` の条件が、**すべての** pre-Frobenioid で自動的に成り立つ

可能性を排除できない。★その場合、条件を書いた意味が無い。

★これは我々自身のゲートの **D12「load-bearing なのに負の対照が無い」** と同じ形の穴である
(`Found/` には D12 が及ばないのでゲートは PASS するが、規律としては埋める)。

## ★作るもの

`𝒟 = Vee`、`Φ = ℕ` はそのままに、`𝒞` を **`𝔽_Φ` の「次数 1 の射」だけからなる**
**広い部分圏**に取る。すると

* `PreFrobenioid` の4条件は**すべて満たす**(下で証明)
* しかし **(ii)(各次数 `n` の Frobenius 型射の存在)が落ちる** ——
  `n = 2` の射がそもそも1本も無い

★対象は減らさず**射だけ**を減らすので、`totally epimorphic` は遺伝する。
これが `WideSubcategory` を使う理由である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

/-! ### 次数 1 の射がなす広い部分圏 -/

/-- 「Frobenius 次数が 1」という射のクラス。 -/
def deg1Prop : MorphismProperty wC := fun _ _ φ => wn φ = 1

instance : deg1Prop.ContainsIdentities := ⟨fun _ => rfl⟩

instance : deg1Prop.IsStableUnderComposition :=
  ⟨fun f g hf hg => by
    show wn (f ≫ g) = 1
    rw [wn_comp, show wn g = 1 from hg, show wn f = 1 from hf, mul_one]⟩

instance : MorphismProperty.IsMultiplicative deg1Prop where

/-- ★負の対照の `𝒞` —— `𝔽_Φ` の次数 1 の射だけを残した広い部分圏。 -/
abbrev nC : Type := WideSubcategory deg1Prop

/-! ### `PreFrobenioid` の4条件を確かめる -/

/-- ★**広い部分圏は totally epimorphic 性を継ぐ**。

包含は忠実なので、`f ≫ g = f ≫ h` を下の圏へ落として `wC` の epi 性を使えばよい。 -/
theorem isTotallyEpimorphic_nC : IsTotallyEpimorphic nC := by
  intro X Y f
  refine ⟨fun {Z} g h hgh => ?_⟩
  refine InducedWideCategory.Hom.ext ?_
  haveI : Epi f.hom := wP.totEpiC _ _ _
  refine (cancel_epi f.hom).mp ?_
  calc f.hom ≫ g.hom = (f ≫ g).hom := rfl
    _ = (f ≫ h).hom := by rw [hgh]
    _ = f.hom ≫ h.hom := rfl

/-- ★**`nC` は connected**(2026-08-16 に追加)。

★次数 1 に絞っても `Vee` の射の持ち上げ `⟨f, 0, 1⟩` は残るので、
連結性は落ちない。★**(ii) が落ちるのは次数 2 が無いからであって、
連結性とは無関係である**ことがここで見える。 -/
instance isConnected_nC : IsConnected nC := by
  refine IsConnected.of_induct (j₀ := (⟨⟨Vee.top⟩⟩ : nC)) ?_
  intro p hp0 hstep A
  have key : ∀ d : Vee, (⟨⟨d⟩⟩ : nC) ∈ p :=
    induct_on_objects (J := Vee) {d | (⟨⟨d⟩⟩ : nC) ∈ p} hp0
      (fun {d₁ d₂} f => hstep (⟨⟨f, 0, 1⟩, rfl⟩ :
        (⟨⟨d₁⟩⟩ : nC) ⟶ (⟨⟨d₂⟩⟩ : nC)))
  exact key A.obj.base

/-- ★**負の対照の pre-Frobenioid**。`Definition 1.1, (iv)` の条件をすべて満たす。

* 関手 `𝒞 → 𝔽_Φ` = 広い部分圏の包含
* `Φ` は divisorial(`isDivisorial_nat`)
* `𝒞` は totally epimorphic(上で証明)
* `𝒟 = Vee` は totally epimorphic(`isTotallyEpimorphic_vee`)
* ★`𝒞`, `𝒟` は connected(2026-08-16 に追加)
-/
def nP : PreFrobenioid nC wΦ where
  toElem := wideSubcategoryInclusion deg1Prop
  divisorial _ := isDivisorial_nat
  totEpiC := isTotallyEpimorphic_nC
  totEpiD := isTotallyEpimorphic_vee
  connectedC := isConnected_nC
  connectedD := isConnected_vee

/-- 包含なので `deg_Fr` は下の圏の次数そのもの。 -/
@[simp] theorem nP_degFr {A B : nC} (φ : A ⟶ B) : nP.degFr φ = wn φ.hom := rfl

/-- ★**`nP` の射はすべて次数 1** —— これが (ii) を落とす。 -/
theorem nP_degFr_eq_one {A B : nC} (φ : A ⟶ B) : nP.degFr φ = 1 := φ.property

/-! ### ★★判定 —— `nP` は Frobenioid ではない -/

/-- ★★**`nP` は Frobenioid ではない。**

落ちるのは **(ii)**(`frobDegSurj`): 次数 2 の Frobenius 型射を要求するが、
`nC` には次数 2 の射が**1本も無い**。

★したがって `Definition 1.3` は `PreFrobenioid` の**真の**制限である ——
条件は自動的には成り立たない。 -/
theorem nNot_frobenioid : ¬ Frobenioid nP := by
  intro F
  obtain ⟨B, φ, -, hdeg⟩ := F.core.frobDegSurj ⟨wTop⟩ 2
  have h1 : nP.degFr φ = 1 := nP_degFr_eq_one φ
  rw [hdeg] at h1
  exact absurd h1 (by decide)

/-- ★**落ちる条を名指しする** —— (ii) の `frobDegSurj` が `nP` で成り立たない。

★`nNot_frobenioid` は「Frobenioid でない」を言うだけだが、こちらは
**どの条が落ちるか**を機械で固定する。 -/
theorem nNot_frobDegSurj :
    ¬ (∀ (A : nC) (n : ℕ+), ∃ (B : nC) (φ : A ⟶ B),
        IsFrobeniusType nP φ ∧ nP.degFr φ = n) := by
  intro h
  obtain ⟨B, φ, -, hdeg⟩ := h ⟨wTop⟩ 2
  have h1 : nP.degFr φ = 1 := nP_degFr_eq_one φ
  rw [hdeg] at h1
  exact absurd h1 (by decide)

/-- ★★**`Definition 1.3` は `PreFrobenioid` の真の制限である。**

満たす例(`wP`)と満たさない例(`nP`)の**両方**が存在する。 -/
theorem frobenioid_is_proper_restriction :
    (∃ (P : PreFrobenioid wC wΦ), Frobenioid P) ∧
      (∃ (P : PreFrobenioid nC wΦ), ¬ Frobenioid P) :=
  ⟨⟨wP, wIsFrobenioid⟩, ⟨nP, nNot_frobenioid⟩⟩

/-! ### ★出典の紐付け(`.src`) -/

def nP.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 20, item := "Definition 1.1, (iv) — 負の対照の pre-Frobenioid",
    sectionId := "frdi-def-1-1-iv" }

def nNot_frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Definition 1.3 — 負の対照((ii) が落ちる)",
    sectionId := "frdi-def-1-3-ii" }

end ABC3.Found.FrdI
