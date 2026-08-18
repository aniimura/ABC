---
name: pf-frobenius-absorbed-by-transition
description: [FrdI] 𝒞^pf の要点は「Frobenius 型射を添字の遷移に吸収できる」こと。𝒞 では偽の命題(mono・右から割れる)が 𝒞^pf では真になる。
metadata:
  type: project
---

`Hom^pf(A,B)` は添字圏(遷移＝両脚に Frobenius 型射を合成する)上の filtered colimit。
★★**遷移の第 2 脚に、割りたい Frobenius 型射 `ψ` そのものを置ける**。
すると `idxTransport` の仕様 `φ ≫ ψ = a ≫ (transported)` が使えて、
`a` は pre-Frobenioid の全射性(`totEpiC`)で **epi** だから、両辺から `a` を消せる。

これが `Found/FrdI/Prop32Perfect.lean` の `homPf_cancel_frobType_idx`。

## ★★`𝒞` では偽・`𝒞^pf` では真になる命題

| 命題 | `𝒞` | `𝒞^pf` |
|---|---|---|
| Frobenius 型射は mono | ★**偽**(単元のぶんずれる。捻れ積 `𝔽_Φ ⋉ G` で `u^n=1` なら `φ∘u=φ`) | ★真 `pfRoot_frobTypeMono` |
| Frobenius 型射で**右から割れる** | ★**偽**(`arbFactor` は Frobenius を**左**に出す) | ★真 `pfRoot_frob_div` |
| どの対象も各次数の Frobenius 型射の**終域** | 偽(`frobDegSurj` は始域が自由) | ★真 `pfRoot_frobDegSurj_cod` |

★この 3 つがそれぞれ「忠実・充満・本質的全射」に対応して
**`𝒞^pf ≃ (𝒞^pf)^pf`**(perfection の冪等性)を与える(`toPfRoot_isEquivalence`)。

## ★右から割る手(3 段)

1. 底同型な射は `𝒞` で「Frobenius 型 ≫ pre-step」に分解する
   (`Definition 1.3, (iv)(a)` の 3 分解 → pull-back 部分が底同型 →
   `Proposition 1.4, (iii)` で同型に潰れる)。次数が割り切れるだけなら
   `exists_frob_factor_of_dvd` で Frobenius 部分の次数を `n` に切る。
2. ★★**`(γ₀, χ, ξ)` を 3 脚添字の遷移として使う** —— `γ₀`(分解の Frobenius 部分)と
   `χ`(割る側の代表元)が**同じ次数 `n`** の Frobenius 型射だから正当。
3. `compRoot_mk3` で `𝒞^pf` へ戻す。

**How to apply:** `𝒞^pf` について「`𝒞` で成り立たないから無理」と判断する前に、
**添字の遷移に吸収できないか**を必ず見る。型の罠は [[widesubcategory-type-trap]]、
射の受け方は [[lean-rebind-morphisms-clean-types]]。
