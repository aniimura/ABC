import ABC3.Found.FrdI.Prop114

/-!
# Gap — [FrdI] §1 で原典の仮定からは出なかった 2 段

★**`Found/FrdI/` で §1 の 15 項目のうち 13 項目を完全に実装した。**
残る 2 項目について、**主張を弱めずに**不足を型に出す(G5)。

★★**どちらも「試していない」でも「証明できなかった」でもない** ——
★**不足しているものを式で特定した**段階である。

## ★分類について

`GapClass` の既定は `modelError`(①)であり、`sourceGap`(③)は最後の手段である。
★**ここでは 2 件とも `missingMath`(②)とする** ——
★**反例となる Frobenioid を実際に構成していない**以上、
「原典の飛躍」と名乗るべきではない。
★**`falsifier` には「何が起きれば ① に落ちるか」を書く。**
-/

namespace ABC3.Gap.FrdI

open CategoryTheory ABC3.Found.FrdI

universe v u w u2 v2

/-! ## ★`Proposition 1.14, (iii)` の `⟸`

原文 (FrdI p.41):
> (iii) Suppose that φ is irreducible. Then φ is a non-pre-step if and only if

★**`⟹` は実装済み**(`prop_1_14_iii_mp`、`N = N_𝒟 + 3` で構成的)。

★★**`⟸`(pre-step ⟹ 条件が偽)には「prime-Frobenius 射が FSM 射である」が要る。**
- 条件文は鎖の各因子だけでなく ★**`ψ` 自身にも FSMI を要求する**
  (原文 l.2136「where α1, …, αn, **ψ** are FSMI-morphisms」)
- 合成の次数は乗法的なので、鎖を伸ばすには次数 > 1 の因子が要り、
  (i) によりそれは prime-Frobenius である
- ★**(b) 型((c) 型)の因子は `Div`(`𝒟` の FSMFF (b))で有界**なので、
  他の伸ばし方は無い

★★**そして `Definition 1.3` で mono を与える条項は `preStepMono` **1 つだけ**であり、
pre-step にしか mono を与えない。**

★**穴の正体は `Found/FrdI/Prop114.lean` の 2 本で式にした**:
- `frobType_cancel_invariants` —— Frobenius 型射で消去すると
  **不変量(`Base`・`degFr`・`Div`)はすべて一致する**
- `mono_of_frobType_of_faithful` —— **`𝒞 → 𝔽_Φ` が忠実なら mono**

★★**したがって穴はちょうど「単元の分」だけ開いている** ——
`Definition 1.3` が与えるのは `faithfulUpToUnits`(単元を除く忠実性、
しかも co-angular pre-step についてのみ)である。

★**破れる機構も式にした**(`not_mono_of_frobNormalized_of_torsion`) ——
`Frobenius-normalized` の等式に `𝒪^×` の `p`-捻れ元を当てると
`𝟙 ≫ ζ = α ≫ ζ` になり、`ζ` は mono でない。
-/

/-- ★**`Proposition 1.14, (iii)` の `⟸` に不足しているもの**。

★原典の語彙で書けば「**prime-Frobenius 射は FSM 射である**」。
ここでは mono の部分だけを出す(fiberwise-surjectivity も別途要る)。 -/
structure Gap_1_14_iii {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
    {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) : Prop where
  /-- 不足: **Frobenius 型射が mono である**こと。

  ★`Definition 1.3` の `preStepMono` は pre-step にしか mono を与えない。 -/
  frobTypeMono : ∀ {A B : C} (ε : A ⟶ B), IsFrobeniusType P ε → Mono ε
  /-- 不足: **Frobenius 型射が fiberwise-surjective である**こと。

  ★`Proposition 1.11, (vii)` が与えるのは
  **co-angular pre-step に沿った持ち上げ**だけである。 -/
  frobTypeFS : ∀ {A B : C} (ε : A ⟶ B), IsFrobeniusType P ε → IsFiberwiseSurjective ε

def Gap_1_14_iii.record : ABC3.Meta.GapRecord :=
  { source :=
      { paper := "FrdI", pdfPage := 41, item := "Proposition 1.14, (iii)",
        sectionId := "frdi-prop-1-14" },
    classification := ABC3.Meta.GapClass.missingMath,
    falsifier :=
      "`Definition 1.3` の条項から「Frobenius 型射は mono」または「𝒞 → 𝔽_Φ は忠実」が" ++
      "導かれることが示されれば ① に落ちる。逆に `𝒪^×` に p-捻れを持つ Frobenioid を" ++
      "構成できれば ③ に上がる(機構は `not_mono_of_frobNormalized_of_torsion` で式にした)。" ++
      "なお使用箇所(FrdI p.63 の Theorem 3.4)では追加の仮定は課されていない。" }

/-! ## ★`Proposition 1.6, (v)` の `base-trivial` / `metrically trivial` の `⟸`

原文 (FrdI p.28):
> sub-quasi-Frobenius-trivial; metrically trivial; base-trivial; perfect; group-

★**`⟹` は実装済み**(`cfp_metricallyTrivial_mp` / `cfp_baseTrivial_mp`)。

★★**`⟸` には「`𝒞` の同型を、底を指定して取り直せる」ことが要る。**
`𝒞' = 𝒞 ×_𝒟 𝒟'` の同型 `(f, g)` は四角形
`A.hom ∘ Base f = G(g) ∘ Dd'.hom` を満たさねばならず、
★**`f` の `Base` が 1 つに指定されてしまう。**
`base-trivial` が与えるのは**ある**同型であって、底を指定した同型ではない。

★★**`Proposition 1.6` の仮定は「FSM 射を FSM 射に送る関手 `G : 𝒟' ⥤ 𝒟`」だけ**で、
`Aut-ample` も `G` の充満性も無い。

★**原文が (vi)(`Aut-ample` 等)を「**if**」(片向き)としか書いていないことに注意** ——
★**著者は (v) と (vi) で向きを書き分けている。**
-/

/-- ★**`Proposition 1.6, (v)` の 2 件の `⟸` に不足しているもの**。

★原典の語彙で書けば「**`𝒟` の自己同型が `𝒞` に持ち上がる**」、すなわち
`Definition 1.2, (iv)` の `Aut-ample`。 -/
structure Gap_1_6_v {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
    {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) : Prop where
  /-- 不足: **全対象が `Aut-ample`** であること。

  ★これがあれば「底を指定した同型」を取り直せる。 -/
  autAmple : ∀ A : C, IsAutAmple P A

def Gap_1_6_v.record : ABC3.Meta.GapRecord :=
  { source :=
      { paper := "FrdI", pdfPage := 28, item := "Proposition 1.6, (v)",
        sectionId := "frdi-prop-1-6" },
    classification := ABC3.Meta.GapClass.missingMath,
    falsifier :=
      "CFP の構造(`Definition 1.3` ＋ `G` が FSM 射を保つこと)だけから" ++
      "「底を指定した同型を取り直せる」ことが導かれれば ① に落ちる。" ++
      "逆に `Aut-ample` でない対象を持つ具体的な CFP を構成して反例にできれば ③ に上がる。" ++
      "原文が (vi) を片向きでしか述べていないことは、著者が向きを意識している証拠である。" }

end ABC3.Gap.FrdI
