import ABC3.Found.FrdI.Prop44Core

/-!
# Gap — [FrdI] `Proposition 4.4, (ii)` の最後の 1 条

★**`Found/FrdI/Prop44Core.lean` で `𝒞^birat` の `FrobenioidCore` 21 条のうち
20 条と、`Definition 1.3, (i)(c)` の圏同値(`birat_plBkEquiv`)を証明した。**
★残るのは **`otriBase`(`Definition 1.3, (iii), (c)` の「全単射は `Base(φ)` にしか依らない」)
の 1 条だけ**である。

## ★原文

原文 (FrdI p.85):
> exercise to check that Cbirat is, in fact, a Frobenioid of group-like type. Moreover, it

★原文はここを **「routine exercise」** と書いて証明を置かない。

## ★★測定 —— 残る 1 条は何と同値か(2026-08-17)

`otriBase` は次の形である:

> `φ`・`φ' : A ⟶ B` が co-angular pre-step で `Base φ = Base φ'` のとき、
> `α ∈ 𝒪^▷(A)`・`β ∈ 𝒪^▷(B)` について
> `φ ≫ β = α ≫ φ` ⟹ `φ' ≫ β = α ≫ φ'`。

★★`𝒞^birat` では co-angular pre-step は**同型**なので、`p := φ⁻¹` と置くと
`β = p ≫ α ≫ φ` で、結論は `u ≫ α = α ≫ u`(`u := φ' ≫ p ∈ 𝒪^▷(A)`)に**同値**になる。
★★★すなわち **`otriBase` ⟺ `𝒪^▷(A^birat)`(= `𝒪^×(A^birat)`)が可換**。

★`Remark 1.3.1` は「`𝒪^▷(A)` は可換」を **(iii)(b),(c) から**導くが、
その (c) が `otriBase` 自身なので、**この経路は循環している**。

## ★★測定 —— どこまでは示せたか

★以下は `Prop44Core.lean` で**実際に証明した**もの:

| 事実 | 宣言 |
|---|---|
| `𝒞^birat` の自己射はすべて co-angular | `birat_isCoAngular_endo` |
| よって `𝒪^▷(A^birat)` の元はすべて**同型**(= 群になる) | (上と `birat_isIso_of_coaPre_birat`) |
| `𝒪^▷(A^birat)` の元は `[c]⁻¹ ≫ [s]`(`c`・`s : E ⟶ A` は `𝒞` の co-angular pre-step で `Base c = Base s`)の形 | `birat_hom_repr` |

★★さらに紙の上では次まで出ている(未実装):
**`𝒞` の `otriBase` を `c` と `s` に当てると `Θ_c = Θ_s` なので、
`u = [c]⁻¹ ≫ [s]` は `𝒪^▷(A)_𝒞` の像すべてと可換**である。
★したがって **`𝒪^▷(A)_𝒞` の像は `𝒪^▷(A^birat)` の中心に入る**。

★★★**足りないのは「`𝒪^▷(A)_𝒞` の像とその逆元が `𝒪^▷(A^birat)` を生成する」こと**
—— これが言えれば可換性は上の中心性から**直ちに**従う。
★モデル(数体上の直線束の Frobenioid)では `𝒪^▷(A)^gp = K^×` で成り立つが、
`Definition 1.3` の公理からこれを出す道は見つかっていない
(`preStepSpan` は `X ⟶ A`・`X ⟶ E` の span を与えるだけで、
必要な `A ⟶ E` 向きの co-angular pre-step は与えない)。

## ★★★2026-08-18 の測定 —— 穴は「birat-Frobenius-normalized でない `𝒞`」だけに絞れた

★`Prop25iii.lean` の `otri_mul_comm` は
**「`A` が Frobenius-normalized ⟹ `𝒪^▷(A)` は可換」**を与える。
上で測ったとおり `otriBase` ⟺ `𝒪^▷(A^birat)` が可換なので、

  ★★★**`𝒞` が birationally Frobenius-normalized 型 ⟹ `𝒞^birat` は 21 条をすべて満たす**

(`Found/FrdI/Prop44Otri.lean` の `birat_frobenioidCore_of_frobNormalized`)。

★したがって**穴が残るのは birationally Frobenius-normalized でない `𝒞` だけ**である。
原文 `Example 4.6` はそういう `𝒞` が実在することを示しているので**穴は空虚ではない**が、
`Theorem 5.2` の model Frobenioid・`Example 6.1`・`6.3` はすべて model 型
(したがって birationally Frobenius-normalized 型)なので、
★★**下流で実際に要る場面では本条は埋まっている**。

## ★★★2026-08-18 の訂正 —— 「見る条を間違えていた」

★上の欄に

> `preStepSpan` は `X ⟶ A`・`X ⟶ E` の span を与えるだけで、
> 必要な `A ⟶ E` 向きの co-angular pre-step は与えない

と書いた。★**これは誤りである。**要るものを (i)(b) の `preStepSpan` から出そうとしていたが、
実際に効くのは **(iii)(d) の `coaPreOverEquiv`**、すなわち

  `(𝒞^coa-pre)_A ≃ Order(Φ(A))^opp`

である。★`Order(Φ(A))` は**前順序圏**なので射は高々 1 本、そして `Φ(A)` は
`δ ≼ δ + ε` により**有向**である。ゆえに `(𝒞^coa-pre)_A` は**余フィルター**で、
`δ_c + δ_s` に対応する対象から `(E,c)`・`(E,s)` の双方へ射が降りる。

★★★**したがって Ore の四角形は一般の Frobenioid で作れる**
(`Found/FrdI/Prop44Ore.lean` の `exists_ore_square_coaPre`・`exists_ore_common`、
**証明済み・`sorry` 無し**)。`Thm52Frob.lean` の `exists_ore_square` は model について
因子を明示して作ったものだったが、**model であることは要らなかった**。

## ★★残りを分解した(2026-08-18)

★★**「壁」ではなく道である。**`ResearchPaper/frdi-decomposition.json` の
`otricomm` チェーン(7 節点・5 層)に割った。`node tools/frdi-newleaves.mjs` が層と葉を印字する。

| 層 | 節点 | 状態 |
|---|---|---|
| 0 | `ore-square`(Ore の四角形) | ★**済**(`exists_ore_square_coaPre`) |
| 0 | `divgp-hom`(`𝒪^▷(A^birat) → Φ^gp(A)` が準同型) | ★**済**(`otriDivGpHom`) |
| 0 | `image-central`(`𝒪^▷(A)_𝒞` の像が中心) | ★済(2026-08-17) |
| 1 | `ker-eq-image`(核 = `𝒪^×(A)` の像。(vi) `faithfulUpToUnits` から) | 未 |
| 2 | `central-ext`(中心拡大 ⟹ 交換子が交代双線形形式) | 未 |
| 3 | ★**`pairing-vanishes`(その形式が消えること)** | ★**未。ここが本体** |
| 4 | `otriBase`(可換性からの帰結) | 半済(`birat_otriBase_of_comm`) |

★★★**穴の形が変わった**——「圏論の条件 1 つ」から
「**1 つの交代双線形形式が消えるか**」になった。反例があるならその形式が非零な例として見える。

## ★影響範囲

★`Proposition 4.4, (ii)` が閉じないので、
`Prop 4.4` 全体の条なし `.src` は**まだ書けない**。
★`(i)`・`(iii)` の `Φ^birat`・`(iv)` の辞書は本条に依存しないので、
それらは独立に進められる。
-/

namespace ABC3.Gap.FrdI

/-- ★★**`Proposition 4.4, (ii)` に足りないもの**。

★分類は **② `missingMath`** —— 主張はモデルで真であり、
原典も「routine exercise」と書く。**足りているのは我々の構成**であって、
原典の飛躍だと名乗るには根拠が弱い。 -/
structure Gap_4_4_ii_otriBase : Prop where
  /-- 不足: `𝒪^▷(A)_𝒞` の像とその逆元が `𝒪^×(A^birat)` を生成すること
  (これがあれば中心性から可換性が出る)。 -/
  otriGpGenerates : True

def Gap_4_4_ii_otriBase.record : ABC3.Meta.GapRecord :=
  { source :=
      { paper := "FrdI", pdfPage := 85, item := "Proposition 4.4, (ii)",
        sectionId := "frdi-prop-4-4" },
    classification := ABC3.Meta.GapClass.missingMath,
    falsifier :=
      "★**これが埋まる条件**: `𝒪^▷(A)^gp -> 𝒪^×(A^birat)` が全射であることを " ++
      "`Definition 1.3` の公理から示すこと。★等価な形は " ++
      "「co-angular pre-step `c : E ⟶ A` が与えられたとき、" ++
      "`m ≫ c` と `m ≫ s` がともに `𝒪^▷(A)` に入るような `m : A ⟶ E` を取れる」。" ++
      "★あるいは `𝒪^×(A^birat)` の可換性を直接示す別経路でもよい。" ++
      "★★**③ に上がる条件**: 上のどちらも公理から出ないことを、" ++
      "`Definition 1.3` を満たすが `𝒪^×(A^birat)` が非可換になる Frobenioid の" ++
      "**反例**を作って示すこと。★このとき原文の「routine exercise」は誤りになる。" ++
      "★我々の測定(2026-08-17)では `𝒪^▷(A)_𝒞` の像が中心に入ることまでは出ており、" ++
      "反例があるなら中心拡大の交換子で見えるはずである。" ++
      "★★**2026-08-18 の更新**: (iii)(d) `coaPreOverEquiv` から " ++
      "**Ore の四角形が一般の 𝒞 で作れる**ことを証明した" ++
      "(`Found/FrdI/Prop44Ore.lean` の `exists_ore_square_coaPre` / `exists_ore_common`)。" ++
      "零因子準同型 `otriDivGpHom` も取得した。★残りは " ++
      "`ResearchPaper/frdi-decomposition.json` の `otricomm` チェーンの " ++
      "`ker-eq-image` → `central-ext` → **`pairing-vanishes`** の 3 段であり、" ++
      "最後の 1 段が本体である。" }

end ABC3.Gap.FrdI
