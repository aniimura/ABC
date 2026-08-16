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

- `Definition 1.1, (iii)`(`𝔽_Φ` の合成則、p.19–20)—— ★**一致**
  (`Div(ψ∘φ) = Base(φ)*(Div ψ) + degFr(ψ)·Div(φ)` の向きも確認)
- `Definition 1.1, (iv)`(pre-Frobenioid、p.20)—— ★**一致**。
  要求は「`Φ` が divisorial」「`𝒞`・`𝒟` が connected かつ totally epimorphic」だけで、
  ★**関手 `𝒞 → 𝔽_Φ` に追加条件は無い。**

- `Definition 1.2`(p.21–23)—— ★★**全項一致**:
  - (i) `linear` / `isometric`(★原文どおり **`Div(φ) = 0`**)/ `metrically equivalent`
  - (ii) `base-isomorphism` / `base-isomorphic` / `pull-back morphism`(ファイバー積の形)/
    `base-equivalent` / `base-identity` / **`𝒪^×(A) ⊆ Aut_𝒞(A)`・`𝒪^▷(A) ⊆ End_𝒞(A)`**
    (base-identity かつ linear な自己射)
  - (iii) `pre-step` / `step` / ★**`co-angular`**(3 分解 `φ = α ◦ β ◦ γ` で
    `α` linear・`β` isometric pre-step・`α` か `γ` が base-isomorphism ⟹ `β` が同型)/
    `LB-invertible` / `Frobenius type` / `prime-Frobenius`
  - (iv) `Frobenius-trivial` / `metrically trivial` / `base-trivial` / `isotropic`

★★★**これで `Prop 1.14, (iii)` に入る定義はすべて照合済みになった** ——
`§0`・`Definition 1.1`・`Definition 1.2`・`Definition 1.3`。★**食い違いは 1 つも出なかった。**

★**残る不確かさは「定義の写し」ではなく「命題の証明」の側**である ——
`Proposition 1.4` / `Proposition 1.5` の我々の証明、および捻れ積が 21 条を満たすという
我々の証明。★**後者は機械検証済みなので、定義が正しければ健全である。**

★★**照合中に出た未解決の疑問**(原文側): `Example 3.6` は `Φ = G`(群)を取り
「`𝒞` は group-like 型の Frobenioid」と述べるが、群は sharp でないので
`Definition 1.1, (i)` の意味で **divisorial ではない**。
`Definition 1.1, (iv)` は `Φ` が divisorial であることを要求している。
★**この食い違いを私は解けていない** —— 私の読み違いか、原文の慣用かは未確定。
★なお我々の反例は `Φ^char`(常に divisorial)を divisor monoid に取るので、
★**どちらの読みでも安全側にある。**

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

/-! ## ★★★2026-08-16 —— 穴を 1 点に絞った

★**当初は「`Definition 1.3` からは `Aut-ample` は出ない」と見立てていた。**
反例の設計を 2 度試して 2 度とも `Definition 1.3` に弾かれたので、
★**「出る」方に切り替えて組み立てたところ、ほぼ全部出た。**

`Found/FrdI/Prop110.lean` の `isAutAmple_of_baseTrivial_of_untwisted` が到達点で、
使った道具は:
- `isFrobeniusTrivial_of_baseTrivial`(`Definition 1.3, (i), (a)` ＋ base-trivial)
- `endo_isCoAngular`(`(iii), (b)` を `𝟙` に当てる ⟹ **自己射はすべて co-angular**)
- `preStep_endo_scale`(`(iv), (a)` ＋ `(ii)` ⟹ **`Base` を保って `Div` を `n` 倍**)
- `coaPre_base_diff`(`(iii), (d)` ⟹ **同じ `Div` なら同型 1 本だけ違う**)
- `preStepSpan`(`(i), (b)`)

★★**残る穴はただ 1 つ**:
```
htwist : ∀ s : A ⟶ A, IsPreStep P s → Φ.map (P.Base s) (P.Div s) = P.Div s
```
すなわち ★**pre-step 自己射の底が、自分の `Div` を動かさない**こと。

## ★★穴の性質(2026-08-16 の解析)

`H := Base(Aut_𝒞(A))`、`c(d) := β_d H ∈ Γ/H`(`β_d` は `Div = d` の
co-angular pre-step 自己射の底)と置くと、上の道具から
- `c(n·d) = c(d)`(`preStep_endo_scale`)
- `c(Φ(β_d)(e) + d) = c(d) + c(e)`(合成)
- `c(0) = 0`

が出る。★`Φ(β_d)` が `d` を固定するなら `c(2d) = 2c(d)` と `c(2d) = c(d)` から
★**`c(d) = 0`、すなわち `Aut-ample`**。

★★**`Φ(β_d)` が位数 `m` で `d` を動かすときは `m·c(d) = 0` までしか出ない**
(巡回和 `Σ_{i<m} Φ(β_d)^i(d)` が固定点になることを使う)。
★したがって ★**反例があるとすれば `Γ/H` に位数 `m` の捻れが要る** ——
`Prop 1.14` のときと**同じ形の要求**である。
## ★★★指数 2・`H` 自明の場合は反例が無い(2026-08-16、紙の上で証明)

★`H = 1`、`Γ = ℤ/2 = {1, b}`、`τ := Φ(b)`(`τ² = id`)とする。
`c(d) = 1` なる `d` があると仮定して矛盾を出す:

1. (A) より `c(2d) = c(d) = 1`、また `H = 1` なので `β_{2d} = β_d = b`
2. (B)[`d`, `e = d`] より `c(τ(d) + d) = c(d) + c(d) = 0`。`S := τ(d) + d` と置く
3. `c(S) = 0` かつ `H = 1` なので `β_S = 1`、したがって `Φ(β_S) = id`
4. (B)[`d' = S`, `e = d`] より `c(d + S) = c(S) + c(d) = 1`。
   ★`d + S = τ(d) + 2d`
5. (B)[`d' = 2d`, `e = d`] より `c(τ(d) + 2d) = c(2d) + c(d) = 0`
6. ★★**4 と 5 が同じ元について 1 と 0 を与える —— 矛盾。**

★★**したがって `c ≡ 0`、すなわち `Aut-ample`。**

★**この議論が一般化しない理由もはっきりしている** —— 3 で
「`β_S ∈ H` ⟹ `Φ(β_S) = id`」を使っており、これは
★**`H` が `Φ(A)` に自明に作用する**ことを要求する。
★したがって反例があるとすれば ★★**`Aut_𝒞(A)` の底が `Φ(A)` に非自明に作用する**
ところにしかない。★**探索範囲がさらに狭まった。**(★未決着)

## ★分類

★**`missingMath`(②)のままとする。** ③ を名乗るには反例が要る。
-/

/-- ★**`Proposition 1.6, (v)` の 2 件の `⟸` に不足しているもの**。

★原典の語彙で書けば「**`𝒟` の自己同型が `𝒞` に持ち上がる**」、すなわち
`Definition 1.2, (iv)` の `Aut-ample`。

★★**2026-08-16: これは `htwist` 1 点に絞られた**
(`isAutAmple_of_baseTrivial_of_untwisted` を見よ)。 -/
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
      "★2026-08-16 に穴を 1 点に絞った: `isAutAmple_of_baseTrivial_of_untwisted` は" ++
      "`base-trivial ⟹ Aut-ample` を仮定 `htwist`(pre-step 自己射の底が自分の Div を" ++
      "動かさない)の下で証明する。したがって残るのは `htwist` だけである。" ++
      "★`htwist` が `Definition 1.3` から導かれれば ① に落ち、`Proposition 1.6, (v)` は" ++
      "完全に実装できる。逆に `Γ/H` に捻れを持つ反例を構成できれば ③ に上がる" ++
      "(解析により、反例があるなら `Φ(β_d)` の作用の位数と同じ捻れが `Γ/H` に要る)。" ++
      "原文が (vi) を片向きでしか述べていないことは、著者が向きを意識している証拠である。" }

end ABC3.Gap.FrdI
