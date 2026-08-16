import ABC3.Found.FrdI.Prop114
import ABC3.Check.FrdI.TwistedFrobenioid

/-!
# Gap — [FrdI] §1 で原典の仮定からは出なかった 2 段

★**`Found/FrdI/` で §1 の 15 項目のうち 13 項目を完全に実装した。**
残る 2 項目について、**主張を弱めずに**不足を型に出す(G5)。

★★**どちらも「試していない」でも「証明できなかった」でもない** ——
★**不足しているものを式で特定した**段階である。

## ★分類について

`GapClass` の既定は `modelError`(①)であり、`sourceGap`(③)は最後の手段である。

★★**2026-08-16 に 1 件目を ③(`sourceGap`)へ上げた** ——
★`Check/FrdI/TwistedFrobenioid.lean` に**反例となる Frobenioid を実際に構成した**
(`cx_isFrobenioid`)。★**2 件目は依然 ②(`missingMath`)**である。
★**`falsifier` には「何が起きれば ①② に落ちるか」を書く。**
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

## ★★★反例を構成した(2026-08-16)

★**単なる直積 `𝔽_Φ × G` では駄目である**。合成が
`(f, g) ≫ (f′, g′) = (f ≫ f′, g · g′)` だと、
`ζ ≫ α^d` の `G` 成分は `g_ζ · g_α^d`、`α ≫ ζ` のそれは `g_α · g_ζ` で、
★**`Frobenius-normalized` が `g_α^{d−1} = 1` を強制してしまう**(`G` が自明になる)。

★★**捻れ積 `𝔽_Φ ⋉ G` にすればよい** —— 合成を
```
(f, g) ≫ (f′, g′) = (f ≫ f′, g^{degFr f′} · g′)
```
と定める(`Base`・`Div`・`degFr` は `G` 成分を無視する)。

★★★**`Check/FrdI/TwistedFrobenioid.lean` で、これが `Definition 1.3` の
全条件を満たすことを証明した**(`twFrobenioidCore` 21 条 ＋ `twIsFrobenioid`)。

★**具体化**: `𝒟 = Vee`、`Φ = ℕ`(定数)、`G = ∏_n ℤ/n`。
★★**`G` にすべての位数の捻れがあるので、次数 > 1 の射は 1 本も mono にならない**
(`cx_not_mono_of_deg_gt_one`)。
★`ℤ/2` だけだと「次数 3 の射を使えばよい」という逃げ道が残るが、
★**この `G` はその逃げ道を塞ぐ。**

## ★★★底圏を一対象にすると、主張そのものが反証できる

★★**`𝒟 = Discrete PUnit`(一対象)、`Φ = ℕ`、`G = ∏_n ℤ/n` を取る**
(`Check/FrdI/TwistedFrobenioid.lean` の `cx2P`)。この圏では:

- ★**FSMI 射はちょうど `Div = 1` の射**である ——
  mono から次数 1(`cx2_fsmi_deg`)、irreducible から `Div = 1`(`cx2_fsmi_div`)
- したがって ★★**鎖の長さは `Div` そのもの**(`cx2_chain_div`)
- `Div = 1` の step `cx2Step` は irreducible な **pre-step** だが、
  `cx2Step ≫ ψ` の `Div` は `1 + 1 = 2` で止まる(`cx2_bounded`)

★**原文が (a) の場合に使う「次数を大きくした prime-Frobenius を後置する」手は、
それらが 1 本も mono でないので使えない。**

★★★**したがって `(iii)` の `⟸`(条件 ⟹ 非 pre-step)は偽である**
(`cx2_refutes_1_14_iii`)。
仮定はすべて満たしている —— `Definition 1.3`(`cx2_isFrobenioid`)、
isotropic 型(`cx2_isotropic`)、`𝒟` が FSMFF 型(`cx2_fsmff`)、
`𝒟` が connected・totally epimorphic、`Φ` が divisorial。

★**注意**: 反証は `BoundedFSMIFactor` という**我々の写し方**に対するものである。
原文の条件文は `αn ◦ … ◦ α1 = ψ ◦ φ`(`α` たちと `ψ` が FSMI)の `n` の有界性で、
★**`ψ` に FSMI を課すのは原文の文言そのもの**(上の p.41 の引用)である。

## ★★★逐語照合の結果(2026-08-16)

★**照合した**(原文と 1 条ずつ):
- `§0` の `fiberwise-surjective`(p.14)・`FSM-morphism`・`irreducible`(p.17)・
  `FSMI`・`FSMFF-type`(p.17–18) —— ★**すべて一致**
- `Definition 1.3` の **21 条すべて** —— ★**一致**(合成の向きも確認した)
- (iii) の条件文(p.41)—— `BoundedFSMIFactor` は忠実
- `Frobenius-normalized`(p.23)—— 一致

★**未照合**: `Definition 1.2` の語彙(`pre-step` 等)、`Definition 1.1`。

## ★★★著者は機構を知っている —— `Example 3.6`(p.70)

★★**原文 `Example 3.6` は `G := ℤ ⊕ (⊕_{p∈Primes} ℤ/pℤ)` を取る** ——
★**我々が使ったのと同じ「あらゆる位数の捻れ」**である。そして原文は:

> every morphism of D is either an isomorphism or a non-monomorphism

と述べ、その根拠を **`cf. the existence of the torsion subgroup ⊕_p ℤ/pℤ ⊆ G`**
と明記する。★★**これは我々の `not_mono_of_frobNormalized_of_torsion` そのものである。**

★**ただし原文はそれを「底圏 `𝒟` が FSM 型になる」ことに使う** ——
FSMI 射が 1 本も無くなる方向である。
★★**Frobenioid `𝒞` の側に同じ捻れを置いたとき `Proposition 1.14, (iii)` の
(a) の議論(次数を大きくした prime-Frobenius を後置する)がどうなるかは、
我々が読んだ範囲の本文では論じられていない。**

★したがって現時点の見立ては
「原文の (iii) の (a) の段は、**書かれていない有限性の仮定**を使っている」である。
★**`Definition 1.2` の照合が済むまで、これ以上は言わない。**
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
    classification := ABC3.Meta.GapClass.sourceGap,
    falsifier :=
      "★2026-08-16 に ② から ③ へ上げた。根拠は `Check/FrdI/TwistedFrobenioid.lean` の" ++
      "`cx2_isFrobenioid`(捻れ積 𝔽_{ℕ on 一対象圏} ⋉ ∏_n ℤ/n が Definition 1.3 の" ++
      "全条件を満たす)と `cx2_refutes_1_14_iii`(その圏で irreducible な pre-step が" ++
      "条件を満たしてしまう)。仮定はすべて充足済み(isotropic 型・𝒟 が FSMFF 型・" ++
      "connected・totally epimorphic・Φ が divisorial)。" ++
      "★これが覆るのは (a) この捻れ積が Definition 1.3 のどれかの条項を実は満たさないと" ++
      "示された場合、または (b) 我々の `BoundedFSMIFactor` の写し方が原文と違う場合" ++
      "(どちらも → ①、我々の側の誤り)である。ψ に FSMI を課すのは原文の文言そのもの" ++
      "(p.41「where α1, . . . , αn, ψ are FSMI-morphisms」)であることは確認した。" ++
      "なお使用箇所(FrdI p.63 の Theorem 3.4)では追加の仮定は課されていない。" }

open ABC3.Check.FrdI in
/-- ★★★**`Definition 1.3` は `Gap_1_14_iii` を含意しない** ——
それを満たさない Frobenioid が実在する。

★`cx_isFrobenioid` がその圏が Frobenioid であることを与え、
`cx_frobType_not_mono` が次数 2 の Frobenius 型射が mono でないことを与える。 -/
theorem not_gap_1_14_iii : ¬ Gap_1_14_iii cxP := by
  intro h
  obtain ⟨A, ζ, hft, -, hnm⟩ := cx_frobType_not_mono 2 (le_refl 2)
  exact hnm (h.frobTypeMono ζ hft)

open ABC3.Check.FrdI in
/-- ★★★**`Proposition 1.14, (iii)` の `⟸` そのものが偽である**。

★★**`Gap_1_14_iii` は「原文の証明の一歩が足りない」という記録だったが、
こちらは「主張が成り立たない」という記録である。**
したがって ★**この項目に `.src`(＝完全実装の主張)を付けることはできない。** -/
theorem prop_1_14_iii_mpr_false :
    ¬ (∀ {X : Cx2C} (φ : X ⟶ X), IsIrreducibleMor φ → BoundedFSMIFactor φ →
        ¬ IsPreStep cx2P φ) := by
  intro h
  obtain ⟨X, φ, hirr, hstep, hbdd⟩ := cx2_refutes_1_14_iii
  exact h φ hirr hbdd hstep

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
